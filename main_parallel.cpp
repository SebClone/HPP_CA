#include "hpp_encryptor.hpp"
#include <cmath>
#include <mpi.h>

using Matrix = std::vector<std::vector<uint8_t>>;

/*
    * Wichtige Entscheidungen:
    * Grid wurde in horizontale Streifen aufgeteilt
    * Ghost Zellen im Randbereich sind ein layer groß (jeweils eine Zeile oben und unten)
    * Jeder Prozess hat eine eigene Wandmaske
*/

constexpr int NUM_ITERATIONS = 1000;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    std::string filename = "message.txt";
    std::vector<uint8_t> fileData = readFileBytes(filename); // Read file as binary data
    std::cout << "Read " << fileData.size() << " bytes from " << filename << std::endl;

    size_t originalSize = 0;
    int grid_size = 0;
    if (rank == 0) 
    {
        fileData = readFileBytes(filename);
        originalSize = fileData.size();
        grid_size = static_cast<int>(std::ceil(std::sqrt(originalSize)));
    }

    MPI_Bcast(&originalSize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&grid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int rows_per_rank = grid_size / nprocs;
    int remainder     = grid_size % nprocs;
    int local_rows    = rows_per_rank + (rank < remainder ? 1 : 0);
    int offset_rows   = rank * rows_per_rank + std::min(rank, remainder);

    size_t paddedSize = static_cast<size_t>(grid_size) * grid_size;
    std::vector<uint8_t> padded;
    if (rank == 0) 
    {
        padded = fileData;
        padded.resize(paddedSize, 0x00);
    }

    std::vector<int> sendcounts(nprocs), displs(nprocs);
    if (rank == 0) 
    {
        for (int r = 0; r < nprocs; ++r) 
        {
            int rrows = rows_per_rank + (r < remainder ? 1 : 0);
            sendcounts[r] = rrows * grid_size;
            displs[r]     = ( (r * rows_per_rank) + std::min(r, remainder) ) * grid_size;
        }
    }

    std::vector<uint8_t> local_core(local_rows * grid_size);
    MPI_Scatterv(
        padded.data(), sendcounts.data(), displs.data(), MPI_BYTE,
        local_core.data(), local_rows * grid_size, MPI_BYTE, 
        0, MPI_COMM_WORLD
    );

    Matrix localGrid(local_rows + 2, std::vector<uint8_t>(grid_size));
    for (int i = 0; i < local_rows; ++i)
        for (int j = 0; j < grid_size; ++j)
            localGrid[i+1][j] = local_core[i*grid_size + j];

    Mask wall_mask = generateRandomWallMask(grid_size, 0.1);
    broadcastMask(wall_mask, MPI_COMM_WORLD);

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
    // 1) Ghost-Zeilenaustausch: oben und unten
    MPI_Sendrecv(&localGrid[1][0], grid_size, MPI_BYTE, (rank-1+nprocs)%nprocs, 0,
                 &localGrid[local_rows+1][0], grid_size, MPI_BYTE, (rank+1)%nprocs, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&localGrid[local_rows][0], grid_size, MPI_BYTE, (rank+1)%nprocs, 1,
                 &localGrid[0][0],           grid_size, MPI_BYTE, (rank-1+nprocs)%nprocs, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // 2) Anwenden der Rules auf den Kernbereich (1..local_rows)
    for (int i = 1; i <= local_rows; ++i)
        for (int j = 0; j < grid_size; ++j)
            applyRules(localGrid, wall_mask, true, i, j);

    
    // Extrahiere Kernbereich aus localGrid
    std::vector<uint8_t> result_core(local_rows * grid_size);
    for (int i = 0; i < local_rows; ++i)
        for (int j = 0; j < grid_size; ++j)
            result_core[i*grid_size + j] = localGrid[i+1][j];

    std::vector<uint8_t> result;
    if (rank == 0) result.resize(paddedSize);

    MPI_Gatherv(
        result_core.data(), local_rows * grid_size, MPI_BYTE,
        result.data(), sendcounts.data(), displs.data(), MPI_BYTE,
        0, MPI_COMM_WORLD
    );

    if (rank == 0) {
        // Trimmen auf originalSize und speichern
        result.resize(originalSize);
        saveAsAsciiText(result, "decrypted_message.txt");
    }

    MPI_Finalize();
    return 0;

}