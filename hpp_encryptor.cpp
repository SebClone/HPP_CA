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

void broadcastMask(const Mask& mask, MPI_Comm comm) {
    int rows = mask.size();
    int cols = rows > 0 ? mask[0].size() : 0;
    MPI_Bcast(&rows, 1, MPI_INT, 0, comm);
    MPI_Bcast(&cols, 1, MPI_INT, 0, comm);
    for (int i = 0; i < rows; ++i) {
        MPI_Bcast(const_cast<uint8_t*>(mask[i].data()), cols, MPI_BYTE, 0, comm);
    }
}

void applyRules(std::vector<std::vector<uint8_t>>& grid, const Mask& wall_mask, bool forward)
{
    int grid_size = static_cast<int>(grid.size());

    if (forward)
    {
        // 1) Collision
        for (int i = 0; i < grid_size; ++i)
            for (int j = 0; j < grid_size; ++j)
                grid[i][j] = collision(grid[i][j], wall_mask[i][j]);

        // 2) Propagation
        std::vector<std::vector<uint8_t>> propagation_grid(
            grid_size, std::vector<uint8_t>(grid_size, 0)
        );
        for (int i = 0; i < grid_size; ++i)
        {
            for (int j = 0; j < grid_size; ++j)
            {
                uint8_t center = grid[i][j];
                int up_i    = (i - 1 + grid_size) % grid_size;
                int down_i  = (i + 1) % grid_size;
                int left_j  = (j - 1 + grid_size) % grid_size;
                int right_j = (j + 1) % grid_size;

                uint8_t temp = center;
                uint8_t &up    = propagation_grid[up_i][j];
                uint8_t &down  = propagation_grid[down_i][j];
                uint8_t &left  = propagation_grid[i][left_j];
                uint8_t &right = propagation_grid[i][right_j];

                propagate(temp, up, down, left, right);
                propagation_grid[i][j] |= temp;  // keep any particles that didn't move
            }
        }
        grid = std::move(propagation_grid);

        // 3) Reflection
        for (int i = 0; i < grid_size; ++i)
            for (int j = 0; j < grid_size; ++j)
                grid[i][j] = reflection(grid[i][j], wall_mask[i][j]);
    }
    else
    {
        // 1) Inverse Reflection
        for (int i = 0; i < grid_size; ++i)
            for (int j = 0; j < grid_size; ++j)
                grid[i][j] = inverse_reflection(grid[i][j], wall_mask[i][j]);

        // 2) Inverse Propagation
        // Keep a copy of the post-reflection state
        std::vector<std::vector<uint8_t>> original = grid;
        std::vector<std::vector<uint8_t>> propagation_grid(
            grid_size, std::vector<uint8_t>(grid_size, 0)
        );
        for (int i = 0; i < grid_size; ++i)
        {
            for (int j = 0; j < grid_size; ++j)
            {
                uint8_t center_original = original[i][j];
                uint8_t &center = propagation_grid[i][j];
                // Preserve walls & other upper bits
                center = center_original & 0b11110000;

                int up_i    = (i - 1 + grid_size) % grid_size;
                int down_i  = (i + 1) % grid_size;
                int left_j  = (j - 1 + grid_size) % grid_size;
                int right_j = (j + 1) % grid_size;

                // Note: parameters swapped to reverse the movement
                uint8_t &down  = original[up_i][j];     // north came from below
                uint8_t &up    = original[down_i][j];   // south came from above
                uint8_t &right = original[i][left_j];   // east came from left
                uint8_t &left  = original[i][right_j];  // west came from right

                inverse_propagate(center, up, down, left, right);
            }
        }
        grid = std::move(propagation_grid);

        // 3) Inverse Collision
        for (int i = 0; i < grid_size; ++i)
            for (int j = 0; j < grid_size; ++j)
                grid[i][j] = inverse_collision(grid[i][j], wall_mask[i][j]);
    }
}


