//Implementiert die Funktionen, die in hpp_encryptor.hpp festgelegt sind
#include "hpp_encryptor.hpp"

#include <fstream>
#include <random>

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
    std::ifstream file(filename, std::ios::binary); // open file in binary mode
    std::vector<uint8_t> data;

    if (!file)
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    // Move to the end to get file size
    file.seekg(0, std::ios::end);
    size_t size = file.tellg();

    // Resize vector to fit file contents
    data.resize(size);

    // Go back to start and read data into vector
    file.seekg(0);
    file.read(reinterpret_cast<char *>(data.data()), size);

    return data;
}


std::vector<std::vector<uint8_t>> reshapeToMatrix(const std::vector<uint8_t> &data, size_t &grid_size)
{
    size_t originalSize = data.size();

    // Calculate the smallest square dimension that can hold all the data
    grid_size = static_cast<size_t>(std::ceil(std::sqrt(originalSize)));
    size_t paddedSize = grid_size * grid_size;

    // Copy and pad the data with 0s if needed
    std::vector<uint8_t> padded = data;
    padded.resize(paddedSize, 0);
    // Fill a 2D matrix with the padded data
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
        file << static_cast<char>(byte); // Convert uint8_t to ASCII character
    }

    file.close();
}


Mask generateRandomWallMask(int grid_size, double wall_ratio, uint32_t seed)
{
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

// ----------------------------------------------
// Forward declarations of functions
// ----------------------------------------------
uint8_t collision(uint8_t current_cell, bool is_wall)
{
    if (is_wall)
        return current_cell;

    uint8_t particles = current_cell & 0b00001111;

    if (particles == 0b00000101)
        return (current_cell & 0b11110000) | 0b00001010;
    else if (particles == 0b00001010)
        return (current_cell & 0b11110000) | 0b00000101;

    return current_cell;
}


void propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right)
{
    // Propagation if its an north partical
    if (center & 0b00001000)
    {
        up |= 0b00001000;      // up cell becomes a north partical
        center &= ~0b00001000; // clear the north partical from center
    }
    // Propagation if its an east partical
    if (center & 0b00000100)
    {
        right |= 0b00000100;   // right cell becomes a east partical
        center &= ~0b00000100; // clear the east partical from center
    }
    // Propagation if its an south partical
    if (center & 0b00000010)
    {
        down |= 0b00000010;    // down cell becomes a south partical
        center &= ~0b00000010; // clear the south partical from center
    }
    // Propagation if its an south partical
    if (center & 0b00000001)
    {
        left |= 0b00000001;    // left cell becomes a west partical
        center &= ~0b00000001; // clear the west partical from center
    }
}


uint8_t reflection(uint8_t current_cell, bool is_wall)
{
    if (!is_wall)
        return current_cell;

    uint8_t particles = current_cell & 0b00001111;
    uint8_t reflected = 0;
    if (particles & 0b00001000)
        reflected |= 0b00000010;
    if (particles & 0b00000010)
        reflected |= 0b00001000;
    if (particles & 0b00000100)
        reflected |= 0b00000001;
    if (particles & 0b00000001)
        reflected |= 0b00000100;

    return (current_cell & 0b11110000) | reflected;
}

// ----------------------------------------------
// Backward declarations of functions
// ----------------------------------------------
uint8_t inverse_reflection(uint8_t current_cell, bool is_wall)
{
    return reflection(current_cell, is_wall);
}


void inverse_propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right)
{
    // North particle came from DOWN cell (i.e., came up)
    if (down & 0b00001000)
    {
        center |= 0b00001000;
        down &= ~0b00001000;
    }

    // South particle came from UP cell (i.e., came down)
    if (up & 0b00000010)
    {
        center |= 0b00000010;
        up &= ~0b00000010;
    }

    // East particle came from LEFT
    if (left & 0b00000100)
    {
        center |= 0b00000100;
        left &= ~0b00000100;
    }

    // West particle came from RIGHT
    if (right & 0b00000001)
    {
        center |= 0b00000001;
        right &= ~0b00000001;
    }
}


uint8_t inverse_collision(uint8_t current_cell, bool is_wall)
{
    return collision(current_cell, is_wall);
}

/*
    *-------------------------------------------------------
    * Functions for MPI communication
    * -------------------------------------------------------
*/

void broadcastMask(Mask &mask, MPI_Comm comm) {
    int rows = mask.size();
    int cols = rows > 0 ? mask[0].size() : 0;
    MPI_Bcast(&rows, 1, MPI_INT, 0, comm);
    MPI_Bcast(&cols, 1, MPI_INT, 0, comm);

    // Auf Nicht-Root: Speicher für mask anlegen
    mask.resize(rows);
    for (int i = 0; i < rows; ++i) {
        mask[i].resize(cols);
    }

    for (int i = 0; i < rows; ++i) {
        MPI_Bcast(mask[i].data(), cols, MPI_BYTE, 0, comm);
    }
}

void copy_n_bytes(const uint8_t* src, std::size_t count, uint8_t* dst) {
    for (std::size_t i = 0; i < count; ++i) {
        dst[i] = src[i];
    }
}

uint8_t applyRules(Matrix& grid,
                   const Mask& wall_mask,
                   bool forward,
                   int row,
                   int col)
{
    int size = static_cast<int>(grid.size());

    // Indices mit Wrap-Around
    auto wrap = [&](int x) {
        return (x + size) % size;
    };
    int up_row    = wrap(row - 1);
    int down_row  = wrap(row + 1);
    int left_col  = wrap(col - 1);
    int right_col = wrap(col + 1);

    if (forward)
    {
        // 1) Kollision
        uint8_t cell = grid[row][col];
        uint8_t val  = collision(cell, wall_mask[row][col]);

        // 2) Propagation
        uint8_t temp = val;
        uint8_t up    = grid[up_row][col];
        uint8_t down  = grid[down_row][col];
        uint8_t left  = grid[row][left_col];
        uint8_t right = grid[row][right_col];
        propagate(temp, up, down, left, right);
        val = temp;

        // 3) Reflection
        val = reflection(val, wall_mask[row][col]);
        return val;
    }
    else
    {
        // 1) Inverse Reflection
        uint8_t orig   = grid[row][col];
        uint8_t center = orig & 0b11110000;

        // 2) Inverse Propagation
        uint8_t up_src    = grid[down_row][col];
        uint8_t down_src  = grid[up_row][col];
        uint8_t left_src  = grid[row][right_col];
        uint8_t right_src = grid[row][left_col];
        inverse_propagate(center, up_src, down_src, left_src, right_src);

        // 3) Inverse Collision
        uint8_t val = inverse_collision(center, wall_mask[row][col]);
        return val;
    }
}


