#include "hpp_encryptor.hpp"
#include <cmath>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>

using Matrix = std::vector<std::vector<uint8_t>>;

constexpr int NUM_ITERATIONS = 1000;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Modus: true = Verschlüsseln, false = Entschlüsseln
    bool doEncrypt = true;

    const char *in_filename = doEncrypt ? "message.txt" : "encrypted_message.txt";
    const char *out_filename = doEncrypt ? "encrypted_message.txt" : "decrypted_message.txt";
    const char *key_filename = "wall_mask.key";

    double t_start = MPI_Wtime();

    std::vector<uint8_t> fileData;
    uint64_t originalSize = 0;
    int grid_size = 0;

    if (rank == 0)
    {
        // Einlesen der passenden Eingabedatei je nach Modus
        fileData = readFileBytes(in_filename);
        if (fileData.empty())
        {
            std::cerr << "ERROR: could not read input file '" << in_filename << "'\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        originalSize = fileData.size();
        grid_size = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(originalSize))));
        std::cout << (doEncrypt ? "Encrypting" : "Decrypting") << " mode selected.\n";
        std::cout << "Read " << originalSize << " bytes from " << in_filename
                  << ", grid size = " << grid_size << std::endl;
    }

    MPI_Bcast(&originalSize, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&grid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int rows_per_rank = grid_size / nprocs;
    int remainder = grid_size % nprocs;
    int local_rows = rows_per_rank + (rank < remainder ? 1 : 0);
    int offset_rows = rank * rows_per_rank + std::min(rank, remainder);

    size_t paddedSize = static_cast<size_t>(grid_size) * grid_size;
    std::vector<uint8_t> padded;
    if (rank == 0)
    {
        padded = fileData;
        padded.resize(paddedSize, 0x00);
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

    Matrix current(local_rows + 2, std::vector<uint8_t>(grid_size));
    Matrix nextGrid = current;

    for (int i = 0; i < local_rows; ++i)
    {
        copy_n_bytes(
            local_core.data() + i * grid_size,
            static_cast<std::size_t>(grid_size),
            current[i + 1].data());
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

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter)
    {
        int prev = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
        int next = (rank + 1 < nprocs) ? rank + 1 : MPI_PROC_NULL;

        MPI_Sendrecv(
            current[1].data(), grid_size, MPI_UINT8_T, prev, 0,
            current[local_rows + 1].data(), grid_size, MPI_UINT8_T, next, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(
            current[local_rows].data(), grid_size, MPI_UINT8_T, next, 1,
            current[0].data(), grid_size, MPI_UINT8_T, prev, 1,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 1; i <= local_rows; ++i)
        {
            for (int j = 0; j < grid_size; ++j)
            {
                nextGrid[i][j] = applyRules(current, wall_mask, doEncrypt, i, j);
            }
        }

        // --------------- HIER FRAMES SAMMELN & SPEICHERN ----------------
        // Wir nehmen den "fertigen" Zustand dieser Iteration (nextGrid),
        // schneiden die Halo-Zeilen ab, gathern nach Rank 0 und speichern.
        std::vector<uint8_t> frame_core(local_rows * grid_size);
        for (int i = 0; i < local_rows; ++i)
        {
            // nur echte Zellen (Zeilen 1..local_rows), keine Halos
            copy_n_bytes(nextGrid[i + 1].data(),
                         static_cast<std::size_t>(grid_size),
                         frame_core.data() + i * grid_size);
        }

        std::vector<uint8_t> frame;
        if (rank == 0)
        {
            frame.resize(paddedSize); // kompletter grid_size x grid_size-Frame
        }

        MPI_Gatherv(frame_core.data(), local_rows * grid_size, MPI_UINT8_T,
                    rank == 0 ? frame.data() : nullptr,
                    rank == 0 ? sendcounts.data() : nullptr,
                    rank == 0 ? displs.data() : nullptr,
                    MPI_UINT8_T,
                    0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            save_frame_bin(frame, iter); // -> frames/frame_000123.bin
        }

        std::swap(current, nextGrid);
    }

    std::vector<uint8_t> result_core(local_rows * grid_size);
    for (int i = 0; i < local_rows; ++i)
    {
        copy_n_bytes(
            current[i + 1].data(),
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
        // Auf die Eingangsgröße kürzen (gleich bei Encrypt u. Decrypt)
        result.resize(originalSize);
        saveAsAsciiText(result, out_filename);
        std::cout << "Saved " << out_filename << " (" << originalSize << " bytes)\n";
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