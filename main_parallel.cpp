#include "hpp_encryptor.hpp"
#include <cmath>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

using Matrix = std::vector<std::vector<uint8_t>>;

constexpr int NUM_ITERATIONS = 1000;
constexpr int TAG_NS = 0;

static constexpr const char *ENC_BIN = "encrypted_full.bin";
static constexpr const char *ENC_META = "encrypted_full.meta";
static constexpr const char *PLAIN_IN = "message.txt";
static constexpr const char *IMG_IN = "image_message.png";        // Bild als message
static constexpr const char *DECRYPT_IMG = "decrypted_image.png"; // Decrypted Bild
static constexpr const char *DECRYPT_OUT = "decrypted_message.txt";

#ifndef DUMP_FRAMES
#define DUMP_FRAMES 1
#endif

// Save a gathered N x N grid as a binary PGM (P5) so Python can read it easily.
// Each byte is treated as an 8-bit grayscale pixel.
static void save_frame_pgm(const std::vector<uint8_t> &buf, int N, int iter)
{
    // Construct filename: frames/pgm/frame_0000.pgm
    char name[128];
    std::snprintf(name, sizeof(name), "frames/pgm/frame_%04d.pgm", iter);
    std::ofstream out(name, std::ios::binary);
    if (!out)
    {
        std::cerr << "ERROR: could not open " << name << " for writing\n";
        return;
    }
    // PGM header
    out << "P5\n"
        << N << " " << N << "\n255\n";
    // If the buffer size is smaller than N*N (shouldn't happen here), pad with zeros.
    if (buf.size() < static_cast<size_t>(N) * static_cast<size_t>(N))
    {
        out.write(reinterpret_cast<const char *>(buf.data()), buf.size());
        std::vector<uint8_t> pad(static_cast<size_t>(N) * static_cast<size_t>(N) - buf.size(), 0);
        out.write(reinterpret_cast<const char *>(pad.data()), pad.size());
    }
    else
    {
        out.write(reinterpret_cast<const char *>(buf.data()), static_cast<std::streamsize>(N) * static_cast<std::streamsize>(N));
    }
    out.close();
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Modus: true = Verschlüsseln, false = Entschlüsseln
    bool doEncrypt = false;

    int mode = doEncrypt ? 1 : 0;
    MPI_Bcast(&mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    doEncrypt = (mode != 0);

    const char *key_filename = "wall_mask.key";

    double t_start = MPI_Wtime();

    std::vector<uint8_t> fileData;
    uint64_t originalSize = 0;
    int grid_size = 0;

    if (rank == 0)
    {
        if (doEncrypt)
        {
            // fileData = readFileBytes(PLAIN_IN); // Klartext laden
            fileData = readFileBytes(IMG_IN); // Bild laden
            if (fileData.empty())
            {
                std::cerr << "ERROR: could not read input file '" << PLAIN_IN << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            originalSize = fileData.size();
            grid_size = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(originalSize))));
            std::cout << "Encrypting mode selected.\n";
            std::cout << "Read " << originalSize << " bytes from " << PLAIN_IN
                      << ", grid size = " << grid_size << std::endl;
        }
        else
        {
            // Meta laden (liefert Originalgröße & N)
            uint32_t Nmeta = 0;
            if (!loadEncryptedMeta(originalSize, Nmeta, ENC_META))
            {
                std::cerr << "ERROR: could not read meta file '" << ENC_META << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            grid_size = static_cast<int>(Nmeta);

            // Voll verschlüsseltes Feld (N*N) laden
            fileData = readFileBytes(ENC_BIN);
            const size_t expected = static_cast<size_t>(grid_size) * grid_size;
            if (fileData.size() != expected)
            {
                std::cerr << "ERROR: '" << ENC_BIN << "' has " << fileData.size()
                          << " bytes, expected " << expected << "\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            std::cout << "Decrypting mode selected.\n";
            std::cout << "Loaded " << fileData.size() << " bytes from " << ENC_BIN
                      << " and meta from " << ENC_META << ", grid size = " << grid_size
                      << ", originalSize = " << originalSize << std::endl;
        }
    }

    MPI_Bcast(&originalSize, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&grid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (nprocs > grid_size)
    {
        if (rank == 0)
            std::cerr << "nprocs > grid_size: mindestens ein Rank bekäme 0 Zeilen – für Torus ungültig.\n";
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    int rows_per_rank = grid_size / nprocs;
    int remainder = grid_size % nprocs;
    int local_rows = rows_per_rank + (rank < remainder ? 1 : 0);
    int offset_rows = rank * rows_per_rank + std::min(rank, remainder);

    size_t paddedSize = static_cast<size_t>(grid_size) * grid_size;
    std::vector<uint8_t> padded;
    if (rank == 0)
    {
        if (doEncrypt)
        {
            padded = fileData;
            padded.resize(paddedSize, 0x00); // Nur beim Encrypt wird gepaddet
        }
        else
        {
            // Beim Decrypt ist fileData bereits exakt N*N groß
            padded = std::move(fileData); // keine Nullen anhängen!
        }
    }

    std::vector<int> sendcounts, displs;
    if (rank == 0)
    {
        sendcounts.resize(nprocs);
        displs.resize(nprocs);
        for (int r = 0; r < nprocs; ++r)
        {
            int rrows = rows_per_rank + (r < remainder ? 1 : 0);
            sendcounts[r] = rrows * grid_size;
            displs[r] = (r * rows_per_rank + std::min(r, remainder)) * grid_size;
        }
    }

    std::vector<uint8_t> local_core(local_rows * grid_size);
    MPI_Scatterv(
        rank == 0 ? padded.data() : nullptr,
        rank == 0 ? sendcounts.data() : nullptr,
        rank == 0 ? displs.data() : nullptr,
        MPI_UINT8_T,
        local_core.data(),
        local_rows * grid_size,
        MPI_UINT8_T,
        0,
        MPI_COMM_WORLD);

    Matrix grid_current(local_rows + 2, std::vector<uint8_t>(grid_size));
    Matrix grid_next = grid_current;

    for (int i = 0; i < local_rows; ++i)
    {
        copy_n_bytes(
            local_core.data() + i * grid_size,
            static_cast<std::size_t>(grid_size),
            grid_current[i + 1].data());
    }

    Mask wall_mask;
    if (rank == 0)
    {
        if (doEncrypt)
        {
            wall_mask = generateRandomWallMask(grid_size, 0.1);
            saveWallMaskBinary(wall_mask, key_filename);
            std::cout << "Wall mask generated and saved to " << key_filename << "\n";
        }
        else
        {
            wall_mask = loadWallMaskBinary(grid_size, key_filename);
            std::cout << "Wall mask loaded from " << key_filename << "\n";
        }
    }
    broadcastMask(wall_mask, MPI_COMM_WORLD);

    // MPI-Topologie für Nachbarschafts-Kommunikation
    int dims[1] = {nprocs}; // 1D-Topologie mit nprocs Prozessen
    int periods[1] = {1};   // 1 = periodisch (Torus), 0 = nicht periodisch
    int reorder = 0;        // MPI darf Ranks nicht neu anordnen für bessere Nachbarschaft
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &cart_comm);

    int up, down;
    MPI_Cart_shift(cart_comm, 0, 1, &up, &down);

    MPI_Request halo_reqs_current[4], halo_reqs_next[4];

    MPI_Recv_init(grid_current[0].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_current[0]);
    MPI_Recv_init(grid_current[local_rows + 1].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_current[1]);
    MPI_Send_init(grid_current[1].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_current[2]);
    MPI_Send_init(grid_current[local_rows].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_current[3]);

    MPI_Recv_init(grid_next[0].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_next[0]);
    MPI_Recv_init(grid_next[local_rows + 1].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_next[1]);
    MPI_Send_init(grid_next[1].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_next[2]);
    MPI_Send_init(grid_next[local_rows].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_next[3]);

    Matrix *active_grid = &grid_current;
    Matrix *target_grid = &grid_next;
    MPI_Request *active_reqs = halo_reqs_current;
    MPI_Request *inactive_reqs = halo_reqs_next;

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter)
    {
        MPI_Startall(4, active_reqs);

        // Innenbereich berechnen
        if (local_rows >= 3)
        {
            for (int i = 2; i <= local_rows - 1; i++)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    (*target_grid)[i][j] = applyRules(*active_grid, wall_mask, doEncrypt, i, j, offset_rows);
                }
            }
        }

        MPI_Waitall(4, active_reqs, MPI_STATUSES_IGNORE);

        if (local_rows >= 1)
        {
            for (int j = 0; j < grid_size; ++j)
            {
                (*target_grid)[1][j] = applyRules(*active_grid, wall_mask, doEncrypt, 1, j, offset_rows);
            }
        }
        if (local_rows >= 2)
        {
            for (int j = 0; j < grid_size; ++j)
            {
                (*target_grid)[local_rows][j] = applyRules(*active_grid, wall_mask, doEncrypt, local_rows, j, offset_rows);
            }
        }
        std::swap(active_grid, target_grid);
        std::swap(active_reqs, inactive_reqs);

// Frames speichern
#if DUMP_FRAMES
        if (iter % 5 == 0 || iter == NUM_ITERATIONS - 1)
        {
            // Gather das aktuelle Grid zu Rank 0
            std::vector<uint8_t> frame_core(local_rows * grid_size);
            for (int i = 0; i < local_rows; ++i)
            {
                copy_n_bytes(
                    (*active_grid)[i + 1].data(),
                    static_cast<std::size_t>(grid_size),
                    frame_core.data() + i * grid_size);
            }

            std::vector<uint8_t> frame;
            if (rank == 0)
            {
                frame.resize(paddedSize);
            }

            MPI_Gatherv(frame_core.data(), local_rows * grid_size, MPI_UINT8_T,
                        rank == 0 ? frame.data() : nullptr,
                        rank == 0 ? sendcounts.data() : nullptr,
                        rank == 0 ? displs.data() : nullptr,
                        MPI_UINT8_T,
                        0, MPI_COMM_WORLD);

            if (rank == 0)
            {
                save_frame_bin(frame, iter);
                // Also save a viewable grayscale frame for visualization/simulation in Python.
                save_frame_pgm(frame, grid_size, iter);
            }
        }
#endif
    }
    for (auto &r : halo_reqs_current)
        MPI_Request_free(&r);
    for (auto &r : halo_reqs_next)
        MPI_Request_free(&r);
    MPI_Comm_free(&cart_comm);

    std::vector<uint8_t> result_core(local_rows * grid_size);
    for (int i = 0; i < local_rows; ++i)
    {
        copy_n_bytes(
            (*active_grid)[i + 1].data(),
            static_cast<std::size_t>(grid_size),
            result_core.data() + i * grid_size);
    }

    std::vector<uint8_t> result;
    if (rank == 0)
    {
        result.resize(paddedSize);
    }
    MPI_Gatherv(
        result_core.data(), local_rows * grid_size, MPI_UINT8_T,
        rank == 0 ? result.data() : nullptr,
        rank == 0 ? sendcounts.data() : nullptr,
        rank == 0 ? displs.data() : nullptr,
        MPI_UINT8_T,
        0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        if (doEncrypt)
        {
            // VOLLEN Zustand + Meta speichern
            saveBinary(result, ENC_BIN);
            saveEncryptedMeta(originalSize, static_cast<uint32_t>(grid_size), ENC_META);
            std::cout << "Saved " << ENC_BIN << " (" << result.size() << " bytes) and "
                      << ENC_META << " (originalSize=" << originalSize
                      << ", grid_size=" << grid_size << ")\n";
        }
        else
        {
            // // Rückgewonnenen Klartext auf Originalgröße kürzen und als Text speichern
            // result.resize(originalSize);
            // saveAsAsciiText(result, DECRYPT_OUT);
            // std::cout << "Saved " << DECRYPT_OUT << " (" << originalSize << " bytes)\n";

            result.resize(originalSize);     // resize to original size
            saveBinary(result, DECRYPT_IMG); // write exact bytes back as binary (PNG/JPG)
            std::cout << "Saved " << DECRYPT_IMG << " (" << originalSize << " bytes)\n";
        }
    }

    double t_end = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << std::fixed << std::setprecision(6)
                  << "Total runtime: " << (t_end - t_start) << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}
