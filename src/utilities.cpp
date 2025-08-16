#include "utilities.hpp"

#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

void printBits(uint8_t value)
{
    std::bitset<8> bits(value);
    std::cout << bits;
}

void printGrid(const std::vector<std::vector<uint8_t>> &grid)
{
    for (const auto &row : grid)
    {
        for (const auto &cell : row)
        {
            printBits(cell);
            std::cout << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<uint8_t> readFileBytes(const std::string &filename)
{
    std::ifstream file(filename, std::ios::binary);
    std::vector<uint8_t> data;

    if (!file)
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    file.seekg(0, std::ios::end);
    size_t size = file.tellg();

    data.resize(size);

    file.seekg(0);
    file.read(reinterpret_cast<char *>(data.data()), size);

    return data;
}

std::vector<std::vector<uint8_t>> reshapeToMatrix(const std::vector<uint8_t> &data, size_t &grid_size)
{
    size_t originalSize = data.size();

    grid_size = static_cast<size_t>(std::ceil(std::sqrt(originalSize)));
    size_t paddedSize = grid_size * grid_size;

    std::vector<uint8_t> padded = data;
    padded.resize(paddedSize, 0);

    std::vector<std::vector<uint8_t>> matrix(grid_size, std::vector<uint8_t>(grid_size));
    for (size_t i = 0; i < paddedSize; ++i)
    {
        matrix[i / grid_size][i % grid_size] = padded[i];
    }

    return matrix;
}

std::vector<uint8_t> flattenMatrix(const std::vector<std::vector<uint8_t>> &matrix)
{
    std::vector<uint8_t> flat;
    for (const auto &row : matrix)
    {
        for (uint8_t value : row)
        {
            flat.push_back(value);
        }
    }
    return flat;
}

void saveAsAsciiText(const std::vector<uint8_t> &data, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (uint8_t byte : data)
    {
        file << static_cast<char>(byte);
    }

    file.close();
}

void save_frame_bin(const std::vector<uint8_t> &frame, int iter)
{
    std::ostringstream oss;
    oss << "frames/frame_" << std::setw(6) << std::setfill('0') << iter << ".bin";

    std::string filename = oss.str();
    std::ofstream ofs(filename, std::ios::binary);

    if (!ofs)
    {
        std::cerr << "ERROR: could not open file '" << filename << "' for writing.\n";
        return;
    }

    ofs.write(reinterpret_cast<const char *>(frame.data()), frame.size());

    if (!ofs)
    {
        std::cerr << "ERROR: failed to write frame to '" << filename << "'\n";
    }
}

Mask generateRandomWallMask(int grid_size, double wall_ratio, uint32_t seed)
{
    if (seed == 0)
    {
        seed = static_cast<uint32_t>(std::random_device{}());
    }
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    Mask wall_mask(grid_size, std::vector<uint8_t>(grid_size, 0));
    int wall_count = 0;

    for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
        {
            if (dist(rng) < wall_ratio)
            {
                wall_mask[i][j] = 1;
                wall_count++;
            }
        }
    }

    std::cout << "Wall cells placed: " << wall_count << std::endl;
    return wall_mask;
}

void saveWallMaskBinary(const Mask &wall_mask, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Cannot write wall mask to " << filename << std::endl;
        return;
    }

    for (const auto &row : wall_mask)
    {
        for (uint8_t wall : row)
        {
            file.write(reinterpret_cast<char *>(&wall), 1);
        }
    }
}

Mask loadWallMaskBinary(int grid_size, const std::string &filename)
{
    Mask wall_mask(grid_size, std::vector<uint8_t>(grid_size, 0));
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Cannot read wall mask from " << filename << std::endl;
        return wall_mask;
    }

    for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
        {
            char wall;
            file.read(&wall, 1);
            wall_mask[i][j] = (wall == 1);
        }
    }

    return wall_mask;
}

void broadcastMask(Mask &mask, MPI_Comm comm)
{
    int rows = mask.size();
    int cols = rows > 0 ? mask[0].size() : 0;
    MPI_Bcast(&rows, 1, MPI_INT, 0, comm);
    MPI_Bcast(&cols, 1, MPI_INT, 0, comm);

    mask.resize(rows);
    for (int i = 0; i < rows; ++i)
    {
        mask[i].resize(cols);
    }

    for (int i = 0; i < rows; ++i)
    {
        MPI_Bcast(mask[i].data(), cols, MPI_BYTE, 0, comm);
    }
}

void copy_n_bytes(const uint8_t *src, std::size_t count, uint8_t *dst)
{
    for (std::size_t i = 0; i < count; ++i)
    {
        dst[i] = src[i];
    }
}

void saveBinary(const std::vector<uint8_t> &data, const char *filename)
{
    std::ofstream out(filename, std::ios::binary);
    out.write(reinterpret_cast<const char *>(data.data()), static_cast<std::streamsize>(data.size()));
}

void saveEncryptedMeta(uint64_t originalSize,
                       uint32_t gridSize,
                       uint64_t startOffset,
                       const char *metaFilename)
{
    std::ofstream out(metaFilename, std::ios::binary);
    const uint32_t magic = 0x48505031; // "HPP1"
    const uint32_t version = 2;        // Version 2: mit StartOffset
    out.write(reinterpret_cast<const char *>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char *>(&version), sizeof(version));
    out.write(reinterpret_cast<const char *>(&originalSize), sizeof(originalSize));
    out.write(reinterpret_cast<const char *>(&gridSize), sizeof(gridSize));
    out.write(reinterpret_cast<const char *>(&startOffset), sizeof(startOffset));
}

bool loadEncryptedMeta(uint64_t &originalSize,
                       uint32_t &gridSize,
                       uint64_t &startOffset,
                       const char *metaFilename)
{
    std::ifstream in(metaFilename, std::ios::binary);
    if (!in)
        return false;
    uint32_t magic = 0, version = 0;
    in.read(reinterpret_cast<char *>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char *>(&version), sizeof(version));
    if (magic != 0x48505031 || version != 1)
        return false;
    in.read(reinterpret_cast<char *>(&originalSize), sizeof(originalSize));
    in.read(reinterpret_cast<char *>(&gridSize), sizeof(gridSize));
    return static_cast<bool>(in);
    if (version >= 2)
    {
        in.read(reinterpret_cast<char *>(&startOffset), sizeof(startOffset));
    }
    else
    {
        startOffset = 0; // Abwärtskompatibilität: Version 1 kennt das Feld nicht
    }
}
