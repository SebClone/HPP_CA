// see folder /include
#include "hpp_rules.hpp"
#include "io.hpp"
#include "utilities.hpp"
#include "config.hpp"

// relevant imports (include what you use principle)
#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <type_traits> // std::true_type // std::false_type

// special datatypes for the grid and the wall mask
using Matrix = std::vector<std::vector<uint8_t>>;
using Mask = std::vector<std::vector<uint8_t>>;

// Specific communication tags for halo transfer (eliminate risk for duplicate meaning)
constexpr int TAG_FROM_UP_A = 100;
constexpr int TAG_FROM_DOWN_A = 101;
constexpr int TAG_FROM_UP_B = 102;
constexpr int TAG_FROM_DOWN_B = 103;

int main(int argc, char **argv)
{
    // MPI setup
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // timers
    double t_io_read = 0.0;
    double t_mask = 0.0;
    double t_loop = 0.0;
    double t_io_write = 0.0;

    // Read config and save to variables
    auto cfg = get_config();
    cfg.mode = parse_mode_from_cli(argc, argv, cfg.mode);
    cfg.grid_size = parse_grid_from_cli(argc, argv, cfg.grid_size == 0 ? 50 : cfg.grid_size);
    cfg.iterations = parse_iters_from_cli(argc, argv, cfg.iterations);

    std::string err;
    if (!validate_config(cfg, err))
    {
        if (rank == 0)
            std::cerr << "ERROR: " << err << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const bool doEncrypt = (cfg.mode == AppMode::Encrypt);
    const int numIterations = cfg.iterations;
    const int frameInterval = cfg.frame_interval;
    const bool dumpFrames = cfg.dump_frames;
    const double wallDensity = cfg.wall_density;
    const std::uint64_t seed = cfg.seed;

    const std::string &inPlain = cfg.input;
    const std::string &encBin = cfg.enc_bin;
    const std::string &metaPath = cfg.meta;
    const std::string &keyPath = cfg.key;
    const std::string &outPlain = cfg.output;

    // Print configuration if needed
    /*
    if (rank == 0)
    {
        std::cout << "Configuration:\n";
        std::cout << "Mode: " << (doEncrypt ? "Encrypt" : "Decrypt") << "\n";
        std::cout << "Iterations: " << numIterations << "\n";
        std::cout << "Grid Size: " << cfg.grid_size << "\n";
        std::cout << "Frame Interval: " << frameInterval << "\n";
        std::cout << "Dump Frames: " << (dumpFrames ? "Yes" : "No") << "\n";
        std::cout << "Wall Density: " << wallDensity * 100.0 << "%\n";
        std::cout << "Seed: " << seed << "\n";
        std::cout << "Input File: " << inPlain << "\n";
        std::cout << "Encrypted Binary File: " << encBin << "\n";
        std::cout << "Meta File: " << metaPath << "\n";
        std::cout << "Key File: " << keyPath << "\n";
        std::cout << "Output File: " << outPlain << "\n";
    }
    */

    double t_start = MPI_Wtime();

    // Initialize metadata
    uint64_t originalSize = 0;
    int grid_size = 0;
    uint64_t start_offset = 0;

    // Initialize MPI datatype for uint64_t (precautionary measure for portability)
    MPI_Datatype MPI_UINT64_MATCHED;
    int rc = MPI_Type_match_size(MPI_TYPECLASS_INTEGER, 8, &MPI_UINT64_MATCHED);
    if (rc != MPI_SUCCESS)
    {
        static_assert(sizeof(unsigned long long) == 8, "Need 64-bit unsigned long long");
        MPI_UINT64_MATCHED = MPI_UNSIGNED_LONG_LONG;
    }

    // Read or initialize metadata (original size, grid size, start offset)
    if (rank == 0)
    {
        uint64_t local_start_offset = 0;
        if (doEncrypt)
        {
            MPI_File fhin;
            if (MPI_File_open(MPI_COMM_SELF, inPlain.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fhin) != MPI_SUCCESS)
            {
                std::cerr << "ERROR: could not open input file '" << inPlain.c_str() << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            MPI_Offset fsz = 0;
            MPI_File_get_size(fhin, &fsz);
            MPI_File_close(&fhin);

            if (fsz < 0)
            {
                std::cerr << "ERROR: invalid size for '" << inPlain.c_str() << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            originalSize = static_cast<uint64_t>(fsz);

            if (cfg.grid_size > 0)
                grid_size = cfg.grid_size;
            else
                grid_size = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(originalSize))));

            const uint64_t cap = static_cast<uint64_t>(grid_size) * grid_size;
            if (originalSize > cap)
            {
                std::cerr << "ERROR: message (" << originalSize << ") exceeds grid capacity (" << grid_size
                          << "x" << grid_size << " = " << cap << ")\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            local_start_offset = (cap - originalSize) / 2;
            std::cout << "Encrypting mode selected.\n";
            std::cout << "Read size " << originalSize << " bytes from " << inPlain.c_str()
                      << ", grid size = " << grid_size << ", start_offset = " << local_start_offset << std::endl;
        }
        else
        {
            uint32_t Nmeta = 0;
            uint64_t meta_start = 0;
            if (!read_meta_rank0(metaPath, originalSize, Nmeta, meta_start))
            {
                std::cerr << "ERROR: could not read meta file '" << metaPath << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            grid_size = static_cast<int>(Nmeta);
            local_start_offset = meta_start;
        }
        start_offset = local_start_offset;
    }

    // Broadcast metadata to all ranks
    MPI_Bcast(&originalSize, 1, MPI_UINT64_MATCHED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&grid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&start_offset, 1, MPI_UINT64_MATCHED, 0, MPI_COMM_WORLD);

    // Abort
    if (nprocs > grid_size)
    {
        if (rank == 0)
            std::cerr << "nprocs > grid_size: at least one rank would get 0 rows – invalid.\n";
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // Data structure for reading/writing in localized 1D chunks (horicontal stripes)
    RowDist dist1D{};
    dist1D.grid_size = grid_size;
    {
        int rows_per_rank = grid_size / nprocs;
        int remainder = grid_size % nprocs;
        dist1D.local_rows = rows_per_rank + (rank < remainder ? 1 : 0);
        dist1D.offset_rows = rank * rows_per_rank + std::min(rank, remainder);
    }
    const std::size_t paddedBytes = static_cast<std::size_t>(grid_size) * grid_size;

    // Rank specific chunk that will be written to or read from
    std::vector<uint8_t> local_core_1d(static_cast<std::size_t>(dist1D.local_rows) * grid_size);

    // Read the initial input (either plain message or binary encrypted data)
    double t_io_r0 = MPI_Wtime();
    if (doEncrypt)
        parallel_read_plain_chunk(inPlain, dist1D, originalSize, start_offset, local_core_1d, false, MPI_COMM_WORLD);
    else
        parallel_read_cipher_chunk(encBin, dist1D, paddedBytes, local_core_1d, false, MPI_COMM_WORLD);
    t_io_read += (MPI_Wtime() - t_io_r0);

    // 2D grid definition, block distribution, 1D -> 2D mapping
    // Define 2D geometry for domain decomposition
    int dims[2] = {0, 0};
    MPI_Dims_create(nprocs, 2, dims);
    int periods[2] = {1, 1}; // Torus topology (wrap-around in up/down and left/right directions)
    int reorder = 0;         
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);

    // Get the neighboring ranks in the 2D geometry
    int up, down, left, right;
    MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
    MPI_Cart_shift(cart_comm, 1, 1, &left, &right);

    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    // 2D block geometry for the grid calculations
    auto split_dim = [](int N, int dim, int coord, int &off, int &len)
    {
        int base = N / dim;
        int rem = N % dim;
        len = base + (coord < rem ? 1 : 0);
        off = coord * base + std::min(coord, rem);
    };
    int off_rows_2d, off_cols_2d, local_rows_2d, local_cols_2d;
    split_dim(grid_size, dims[0], coords[0], off_rows_2d, local_rows_2d);
    split_dim(grid_size, dims[1], coords[1], off_cols_2d, local_cols_2d);

    struct BlockInfo
    {
        int off_r, rows, off_c, cols;
    };
    std::vector<BlockInfo> allBlocks(nprocs);
    std::vector<int> row_off_1d(nprocs), row_len_1d(nprocs);

    for (int r0 = 0; r0 < dims[0]; ++r0)
    {
        int offr, lenr;
        split_dim(grid_size, dims[0], r0, offr, lenr);
        for (int c0 = 0; c0 < dims[1]; ++c0)
        {
            int offc, lenc;
            split_dim(grid_size, dims[1], c0, offc, lenc);
            int cc[2] = {r0, c0};
            int pr;
            MPI_Cart_rank(cart_comm, cc, &pr);
            allBlocks[pr] = {offr, lenr, offc, lenc};
        }
    }
    for (int p = 0; p < nprocs; ++p)
    {
        int r_pr = p;
        int rows_per = grid_size / nprocs;
        int rem = grid_size % nprocs;
        row_len_1d[p] = rows_per + (r_pr < rem ? 1 : 0);
        row_off_1d[p] = r_pr * rows_per + std::min(r_pr, rem);
    }

    // Prepare the grid buffers for the 2D block
    std::vector<uint8_t> gridBufA(static_cast<std::size_t>(local_rows_2d + 2) * (local_cols_2d + 2));
    std::vector<uint8_t> gridBufB(static_cast<std::size_t>(local_rows_2d + 2) * (local_cols_2d + 2));
    auto idx = [W = (local_cols_2d + 2)](int i, int j) -> std::size_t
    {
        return static_cast<std::size_t>(i) * W + static_cast<std::size_t>(j);
    };

    // Send sizes (in bytes) per target rank
    std::vector<int> sendcounts(nprocs, 0), recvcounts(nprocs, 0);
    for (int p = 0; p < nprocs; ++p)
    {
        const auto &B = allBlocks[p];
        int a0 = dist1D.offset_rows, a1 = dist1D.offset_rows + dist1D.local_rows;
        int b0 = B.off_r, b1 = B.off_r + B.rows;
        int s = std::max(a0, b0);
        int e = std::min(a1, b1);
        int rows_overlap = std::max(0, e - s);
        if (rows_overlap > 0)
            sendcounts[p] = rows_overlap * B.cols; // bytes (1 byte per cell)
    }

    // Exchange of sizes
    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Displacements
    auto prefix_sum = [](const std::vector<int> &v)
    {
        std::vector<int> d(v.size(), 0);
        int acc = 0;
        for (size_t i = 0; i < v.size(); ++i)
        {
            d[i] = acc;
            acc += v[i];
        }
        return d;
    };
    std::vector<int> sdispls = prefix_sum(sendcounts);
    std::vector<int> rdispls = prefix_sum(recvcounts);

    int send_total = sdispls.empty() ? 0 : sdispls.back() + sendcounts.back();
    int recv_total = rdispls.empty() ? 0 : rdispls.back() + recvcounts.back();

    std::vector<uint8_t> sendbuf_row2blk(send_total);
    std::vector<uint8_t> recvbuf_row2blk(recv_total);

    // Pack: for each target rank p, all overlapping global rows in ascending order
    for (int p = 0; p < nprocs; ++p)
    {
        const auto &B = allBlocks[p];
        int a0 = dist1D.offset_rows, a1 = dist1D.offset_rows + dist1D.local_rows;
        int b0 = B.off_r, b1 = B.off_r + B.rows;
        int s = std::max(a0, b0);
        int e = std::min(a1, b1);
        int rows_overlap = std::max(0, e - s);
        if (rows_overlap == 0)
            continue;

        int pos = sdispls[p];
        for (int gr = s; gr < e; ++gr)
        {
            std::size_t src_off = static_cast<std::size_t>(gr - dist1D.offset_rows) * grid_size + B.off_c;
            std::memcpy(sendbuf_row2blk.data() + pos,
                        local_core_1d.data() + src_off,
                        static_cast<std::size_t>(B.cols));
            pos += B.cols;
        }
    }

    // Alltoallv
    MPI_Alltoallv(sendbuf_row2blk.data(), sendcounts.data(), sdispls.data(), MPI_BYTE,
                  recvbuf_row2blk.data(), recvcounts.data(), rdispls.data(), MPI_BYTE,
                  MPI_COMM_WORLD);

    // Unpack in gridBufA (without halos, i.e., j-internal starting from 1)
    for (int p = 0; p < nprocs; ++p)
    {
        // Source p owns 1D stripe [row_off_1d[p], row_off_1d[p]+row_len_1d[p])
        int a0 = row_off_1d[p], a1 = row_off_1d[p] + row_len_1d[p];
        int b0 = off_rows_2d, b1 = off_rows_2d + local_rows_2d;
        int s = std::max(a0, b0);
        int e = std::min(a1, b1);
        int rows_overlap = std::max(0, e - s);
        if (rows_overlap == 0)
            continue;

        int pos = rdispls[p];
        for (int gr = s; gr < e; ++gr)
        {
            // local i in the 2D block
            int i_local = (gr - off_rows_2d) + 1;
            std::memcpy(gridBufA.data() + idx(i_local, 1),
                        recvbuf_row2blk.data() + pos,
                        static_cast<std::size_t>(local_cols_2d));
            pos += local_cols_2d;
        }
    }

    // Generate new wall_mask or load existing one
    double t_mask0 = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    Mask wall_mask;
    if (rank == 0)
    {
        if (doEncrypt)
        {
            wall_mask = generateRandomWallMask(grid_size, wallDensity, seed);
            saveWallMaskBinary(wall_mask, keyPath.c_str());
            std::cout << "Wall mask generated and saved to " << keyPath << "\n";
        }
        else
        {
            wall_mask = loadWallMaskBinary(grid_size, keyPath.c_str());
            std::cout << "Wall mask loaded from " << keyPath << "\n";
        }
    }

    // Broadcast the wall mask to all ranks
    std::vector<uint8_t> wall_flat(static_cast<std::size_t>(grid_size) * grid_size);
    if (rank == 0)
    {
        for (int r = 0; r < grid_size; ++r)
            std::memcpy(wall_flat.data() + static_cast<std::size_t>(r) * grid_size,
                        wall_mask[r].data(),
                        static_cast<std::size_t>(grid_size));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(wall_flat.data(), static_cast<int>(wall_flat.size()), MPI_BYTE, 0, MPI_COMM_WORLD);
    if (rank != 0)
    {
        wall_mask.assign(static_cast<std::size_t>(grid_size), std::vector<uint8_t>(grid_size));
        for (int r = 0; r < grid_size; ++r)
            std::memcpy(wall_mask[r].data(),
                        wall_flat.data() + static_cast<std::size_t>(r) * grid_size,
                        static_cast<std::size_t>(grid_size));
    }
    t_mask += (MPI_Wtime() - t_mask0);

    // Define column type for halo exchange
    MPI_Datatype COL_TYPE;
    int rowStride = local_cols_2d + 2; // local columns + 2 for halos
    MPI_Type_vector(local_rows_2d, 1, rowStride, MPI_BYTE, &COL_TYPE);
    MPI_Type_commit(&COL_TYPE);

    // Prepare the grid buffers for the main loop
    uint8_t *active_ptr = gridBufA.data();
    uint8_t *target_ptr = gridBufB.data();

    // ---------------- Main loop ----------------
    // Find out if we are in encryption or decryption mode at compile time
    auto run_main_loop = [&](auto ENC_TAG)
    {
        constexpr bool ENCRYPT = decltype(ENC_TAG)::value;

        if (rank == 0)
        {
            std::cout << "\n[Timing per iteration] (MAX across ranks, seconds):\n";
            std::cout << "iter  comm      inner     border    swap      total\n";
        }

        for (int iter = 0; iter < numIterations; ++iter)
        {
            double it0 = MPI_Wtime();
            double t_comm = 0.0, t_inner = 0.0, t_border = 0.0, t_swap = 0.0;

            // 1) Halo transfers
            double t_post0 = MPI_Wtime();
            MPI_Request reqs[8];
            const bool useA = (active_ptr == gridBufA.data());
            const int TAG_UP = useA ? TAG_FROM_UP_A : TAG_FROM_UP_B;
            const int TAG_DOWN = useA ? TAG_FROM_DOWN_A : TAG_FROM_DOWN_B;
            const int TAG_LEFT = useA ? TAG_FROM_UP_A + 10 : TAG_FROM_UP_B + 10;
            const int TAG_RIGHT = useA ? TAG_FROM_DOWN_A + 10 : TAG_FROM_DOWN_B + 10;

            // RECV halos
            MPI_Irecv(active_ptr + idx(0, 1), local_cols_2d, MPI_BYTE, up, TAG_UP, cart_comm, &reqs[0]);                     // top
            MPI_Irecv(active_ptr + idx(local_rows_2d + 1, 1), local_cols_2d, MPI_BYTE, down, TAG_DOWN, cart_comm, &reqs[1]); // bottom
            MPI_Irecv(active_ptr + idx(1, 0), 1, COL_TYPE, left, TAG_LEFT, cart_comm, &reqs[2]);                             // left
            MPI_Irecv(active_ptr + idx(1, local_cols_2d + 1), 1, COL_TYPE, right, TAG_RIGHT, cart_comm, &reqs[3]);           // right

            // SEND borders
            MPI_Isend(active_ptr + idx(1, 1), local_cols_2d, MPI_BYTE, up, TAG_DOWN, cart_comm, &reqs[4]);             // top inner
            MPI_Isend(active_ptr + idx(local_rows_2d, 1), local_cols_2d, MPI_BYTE, down, TAG_UP, cart_comm, &reqs[5]); // bottom inner
            MPI_Isend(active_ptr + idx(1, 1), 1, COL_TYPE, left, TAG_RIGHT, cart_comm, &reqs[6]);                      // left inner
            MPI_Isend(active_ptr + idx(1, local_cols_2d), 1, COL_TYPE, right, TAG_LEFT, cart_comm, &reqs[7]);          // right inner
            t_comm += (MPI_Wtime() - t_post0);

            // 2) Interior area (without halos)
            double t_in0 = MPI_Wtime();
            if (local_rows_2d >= 3 && local_cols_2d >= 3)
            {
                const int N = grid_size;    // mask width (global)
                const int W = rowStride;    // local stride
                const int j0 = off_cols_2d; // global starting column of our block
#pragma omp parallel for schedule(static)
                for (int i = 2; i <= local_rows_2d - 1; ++i)
                {
                    const int gr = (off_rows_2d + (i - 1)) % N;
                    const int gr_u = (gr - 1 + N) % N;
                    const int gr_d = (gr + 1) % N;

                    // Mask rows: not shifted
                    const uint8_t *__restrict wrow = wall_mask[gr].data();
                    const uint8_t *__restrict wrow_u = wall_mask[gr_u].data();
                    const uint8_t *__restrict wrow_d = wall_mask[gr_d].data();

                    uint8_t *__restrict tgt_row = target_ptr + idx(i, 1);

#pragma omp simd
                    for (int j = 2; j <= local_cols_2d - 1; ++j)
                    {
                        tgt_row[j - 1] = applyRules_fast<ENCRYPT>(
                            active_ptr, W, i, j, wrow, wrow_u, wrow_d, N, j0);
                    }
                }
            }
            t_inner += (MPI_Wtime() - t_in0);

            // 3) Wait + borders
            double t_wait0 = MPI_Wtime();
            MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
            t_comm += (MPI_Wtime() - t_wait0);

            double t_border0 = MPI_Wtime();
            {
                const int N = grid_size;
                const int W = rowStride;
                const int j0 = off_cols_2d;

                // Top row i==1
                if (local_rows_2d >= 1)
                {
                    const int gr = (off_rows_2d + 0) % N;
                    const int gr_u = (gr - 1 + N) % N;
                    const int gr_d = (gr + 1) % N;
                    const uint8_t *__restrict wrow = wall_mask[gr].data();
                    const uint8_t *__restrict wrow_u = wall_mask[gr_u].data();
                    const uint8_t *__restrict wrow_d = wall_mask[gr_d].data();
                    uint8_t *__restrict tgt = target_ptr + idx(1, 1);
#pragma omp simd
                    for (int j = 1; j <= local_cols_2d; ++j)
                    {
                        tgt[j - 1] = applyRules_fast<ENCRYPT>(active_ptr, W, 1, j, wrow, wrow_u, wrow_d, N, j0);
                    }
                }

                // Bottom row i==local_rows_2d
                if (local_rows_2d >= 2)
                {
                    const int iL = local_rows_2d;
                    const int gr = (off_rows_2d + (iL - 1)) % N;
                    const int gr_u = (gr - 1 + N) % N;
                    const int gr_d = (gr + 1) % N;
                    const uint8_t *__restrict wrow = wall_mask[gr].data();
                    const uint8_t *__restrict wrow_u = wall_mask[gr_u].data();
                    const uint8_t *__restrict wrow_d = wall_mask[gr_d].data();
                    uint8_t *__restrict tgt = target_ptr + idx(iL, 1);
#pragma omp simd
                    for (int j = 1; j <= local_cols_2d; ++j)
                    {
                        tgt[j - 1] = applyRules_fast<ENCRYPT>(active_ptr, W, iL, j, wrow, wrow_u, wrow_d, N, j0);
                    }
                }

                // Left column j==1 (excluding corners)
                if (local_cols_2d >= 1 && local_rows_2d >= 3)
                {
                    for (int i = 2; i <= local_rows_2d - 1; ++i)
                    {
                        const int gr = (off_rows_2d + (i - 1)) % N;
                        const int gr_u = (gr - 1 + N) % N;
                        const int gr_d = (gr + 1) % N;
                        const uint8_t *__restrict wrow = wall_mask[gr].data();
                        const uint8_t *__restrict wrow_u = wall_mask[gr_u].data();
                        const uint8_t *__restrict wrow_d = wall_mask[gr_d].data();
                        target_ptr[idx(i, 1)] = applyRules_fast<ENCRYPT>(active_ptr, W, i, 1, wrow, wrow_u, wrow_d, N, j0);
                    }
                }

                // Right column j==local_cols_2d (excluding corners)
                if (local_cols_2d >= 2 && local_rows_2d >= 3)
                {
                    const int jR = local_cols_2d;
                    for (int i = 2; i <= local_rows_2d - 1; ++i)
                    {
                        const int gr = (off_rows_2d + (i - 1)) % N;
                        const int gr_u = (gr - 1 + N) % N;
                        const int gr_d = (gr + 1) % N;
                        const uint8_t *__restrict wrow = wall_mask[gr].data();
                        const uint8_t *__restrict wrow_u = wall_mask[gr_u].data();
                        const uint8_t *__restrict wrow_d = wall_mask[gr_d].data();
                        target_ptr[idx(i, jR)] = applyRules_fast<ENCRYPT>(active_ptr, W, i, jR, wrow, wrow_u, wrow_d, N, j0);
                    }
                }
            }
            t_border += (MPI_Wtime() - t_border0);

            // 4) Swap
            double t_swap0 = MPI_Wtime();
            std::swap(active_ptr, target_ptr);
            t_swap += (MPI_Wtime() - t_swap0);

            double t_iter = MPI_Wtime() - it0;
            t_loop += t_iter;

            // 5) Output (MAX)
            double local_it[5] = {t_comm, t_inner, t_border, t_swap, t_iter};
            double global_it[5] = {0, 0, 0, 0, 0};
            MPI_Reduce(local_it, global_it, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                std::cout << std::fixed << std::setprecision(6)
                          << std::setw(4) << iter << "  "
                          << std::setw(8) << global_it[0] << "  "
                          << std::setw(8) << global_it[1] << "  "
                          << std::setw(8) << global_it[2] << "  "
                          << std::setw(8) << global_it[3] << "  "
                          << std::setw(8) << global_it[4] << "\n";
            }
        }
    };

    // Compile-time dispatch (decision ONCE)
    if (doEncrypt)
        run_main_loop(std::true_type{});
    else
        run_main_loop(std::false_type{});

    // Release datatype
    MPI_Type_free(&COL_TYPE);

    // -------------------------
    // 8) Block(2D) → Row(1D): Alltoallv back, to keep I/O
    // -------------------------
    // result_core_1d: target buffer (without halos), size = dist1D.local_rows * grid_size
    std::vector<uint8_t> result_core_1d(static_cast<std::size_t>(dist1D.local_rows) * grid_size);

    // Send sizes Block->Row (bytes)
    std::vector<int> sendcounts2(nprocs, 0), recvcounts2(nprocs, 0);
    for (int p = 0; p < nprocs; ++p)
    {
        // Target p has row stripe [row_off_1d[p], row_off_1d[p]+row_len_1d[p])
        int a0 = off_rows_2d, a1 = off_rows_2d + local_rows_2d;
        int b0 = row_off_1d[p], b1 = row_off_1d[p] + row_len_1d[p];
        int s = std::max(a0, b0);
        int e = std::min(a1, b1);
        int rows_overlap = std::max(0, e - s);
        if (rows_overlap > 0)
            sendcounts2[p] = rows_overlap * local_cols_2d;
    }
    MPI_Alltoall(sendcounts2.data(), 1, MPI_INT, recvcounts2.data(), 1, MPI_INT, MPI_COMM_WORLD);

    auto prefix_sum2 = [](const std::vector<int> &v)
    {
        std::vector<int> d(v.size(), 0);
        int acc = 0;
        for (size_t i = 0; i < v.size(); ++i)
        {
            d[i] = acc;
            acc += v[i];
        }
        return d;
    };
    std::vector<int> sdispls2 = prefix_sum2(sendcounts2);
    std::vector<int> rdispls2 = prefix_sum2(recvcounts2);
    int send_total2 = sdispls2.empty() ? 0 : sdispls2.back() + sendcounts2.back();
    int recv_total2 = rdispls2.empty() ? 0 : rdispls2.back() + recvcounts2.back();
    std::vector<uint8_t> sendbuf_blk2row(send_total2);
    std::vector<uint8_t> recvbuf_blk2row(recv_total2);

    // Pack: from active_ptr (current buffer after loop) without halos (j=1..local_cols_2d)
    for (int p = 0; p < nprocs; ++p)
    {
        int a0 = off_rows_2d, a1 = off_rows_2d + local_rows_2d;
        int b0 = row_off_1d[p], b1 = row_off_1d[p] + row_len_1d[p];
        int s = std::max(a0, b0);
        int e = std::min(a1, b1);
        int rows_overlap = std::max(0, e - s);
        if (rows_overlap == 0)
            continue;

        int pos = sdispls2[p];
        for (int gr = s; gr < e; ++gr)
        {
            int i_local = (gr - off_rows_2d) + 1;
            std::memcpy(sendbuf_blk2row.data() + pos,
                        active_ptr + idx(i_local, 1),
                        static_cast<std::size_t>(local_cols_2d));
            pos += local_cols_2d;
        }
    }

    MPI_Alltoallv(sendbuf_blk2row.data(), sendcounts2.data(), sdispls2.data(), MPI_BYTE,
                  recvbuf_blk2row.data(), recvcounts2.data(), rdispls2.data(), MPI_BYTE,
                  MPI_COMM_WORLD);

    // Unpack: into result_core_1d at the correct columns (off_cols of the source block)
    for (int p = 0; p < nprocs; ++p)
    {
        // Source p has block allBlocks[p] with column range [off_c .. off_c+cols)
        const auto &SB = allBlocks[p];
        int a0 = row_off_1d[rank], a1 = row_off_1d[rank] + row_len_1d[rank]; // our stripe
        int b0 = SB.off_r, b1 = SB.off_r + SB.rows;                          // their block rows
        int s = std::max(a0, b0);
        int e = std::min(a1, b1);
        int rows_overlap = std::max(0, e - s);
        if (rows_overlap == 0)
            continue;

        int pos = rdispls2[p];
        for (int gr = s; gr < e; ++gr)
        {
            std::size_t dst_off = static_cast<std::size_t>(gr - row_off_1d[rank]) * grid_size + SB.off_c;
            std::memcpy(result_core_1d.data() + dst_off,
                        recvbuf_blk2row.data() + pos,
                        static_cast<std::size_t>(SB.cols));
            pos += SB.cols;
        }
    }

    // -------------------------
    // 9) Write (1D)
    // -------------------------
    double t_io_w0 = MPI_Wtime();
    if (doEncrypt)
    {
        parallel_write_cipher_chunk(encBin, dist1D, result_core_1d, paddedBytes, false, MPI_COMM_WORLD);
        if (rank == 0)
            write_meta_rank0(metaPath, originalSize, static_cast<uint32_t>(grid_size), start_offset);
    }
    else
    {
        parallel_write_plain_trimmed(outPlain, dist1D, result_core_1d, originalSize, start_offset, false, MPI_COMM_WORLD);
    }
    t_io_write += (MPI_Wtime() - t_io_w0);

    // -------------------------
    // 10) Timing Summary
    // -------------------------
    {
        double local_sum[4] = {t_io_read, t_mask, t_loop, t_io_write};
        double global_sum[4] = {0, 0, 0, 0};
        MPI_Reduce(local_sum, global_sum, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "\n[Timing Summary] (MAX across ranks):\n";
            std::cout << "I/O read        : " << global_sum[0] << " s\n";
            std::cout << "Mask+Broadcast  : " << global_sum[1] << " s\n";
            std::cout << "Main loop total : " << global_sum[2] << " s"
                      << "  (avg/iter ~ " << (global_sum[2] / std::max(1, numIterations)) << " s)\n";
            std::cout << "I/O write       : " << global_sum[3] << " s\n";
        }
    }

    double t_end = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << std::fixed << std::setprecision(6)
                  << "Total runtime: " << (t_end - t_start) << " seconds\n";
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
} // main
