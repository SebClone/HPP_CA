#include <iostream>
#include <cstdint> // für uint8_t
#include <bitset>  // für std::bitset
#include <cstring> // for memcpy
#include <fstream>
#include <vector>
#include <random>

const int num_itterations = 10; // Number of iterations for the simulation

// ----------------------------------------------
// Custom functions for printing bits and grid
// ----------------------------------------------
// Prints the 8-bit binary representation of a uint8_t value to std::cout.
// Used for debugging to visualize the bit pattern of a cell.
// Prints the 8-bit binary representation of a uint8_t value to std::cout.
// Used for debugging to visualize the bit pattern of a cell.
// @param value The 8-bit value to print as a binary string.
void printBits(uint8_t value)
{
    std::bitset<8> bits(value);
    std::cout << bits;
}

// Prints the entire 2D grid, showing the bit pattern for each cell.
// Each cell is displayed as 8 bits, separated by spaces.
// Used for debugging and visualization.
// Prints the entire 2D grid, showing the bit pattern for each cell.
// Each cell is displayed as 8 bits, separated by spaces.
// Used for debugging and visualization.
// @param grid The 2D grid to print.
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

// Reads a file as an binary and returns its contents as a vector of bytes
// Reads the contents of a file as binary data and returns it as a vector of bytes.
// @param filename Name of the file to read.
// @return Vector containing the bytes read from the file.
// Reads the contents of a file as binary data and returns it as a vector of bytes.
// @param filename Name of the file to read.
// @return Vector containing the bytes read from the file.
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

// Reshapes a 1D vector into a 2D grid with padding as needed
// Reshapes a 1D vector of bytes into a square 2D grid.
// Pads with zeros if the data does not fit perfectly into a square.
// @param data The input 1D vector of bytes.
// @param grid_size (output) The resulting size of the square grid (number of rows/cols).
// @return 2D vector (grid_size x grid_size) containing the data.
// Reshapes a 1D vector of bytes into a square 2D grid.
// Pads with zeros if the data does not fit perfectly into a square.
// @param data The input 1D vector of bytes.
// @param grid_size (output) The resulting size of the square grid (number of rows/cols).
// @return 2D vector (grid_size x grid_size) containing the data.
std::vector<std::vector<uint8_t>> reshapeToMatrix(const std::vector<uint8_t> &data, int &grid_size)
{
    size_t originalSize = data.size();

    // Calculate the smallest square dimension that can hold all the data
    grid_size = static_cast<int>(std::ceil(std::sqrt(originalSize)));
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
// ----------------------------------------------

// Flattens a 2D matrix back into a 1D vector
// Flattens a 2D grid (matrix) into a 1D vector of bytes in row-major order.
// @param matrix The 2D grid to flatten.
// @return Vector containing the flattened data.
// Flattens a 2D grid (matrix) into a 1D vector of bytes in row-major order.
// @param matrix The 2D grid to flatten.
// @return Vector containing the flattened data.
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

// Saves a vector of bytes as an ASCII text file
// Saves a vector of bytes as an ASCII text file.
// Each byte is written as a character.
// @param data The data to write.
// @param filename The output file name.
// Saves a vector of bytes as an ASCII text file.
// Each byte is written as a character.
// @param data The data to write.
// @param filename The output file name.
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

// ----------------------------------------------
// Wall Bits Mask Functions
// ----------------------------------------------

// Generates a random wall mask for the grid.
// Each cell has a probability `wall_ratio` of being a wall.
// @param grid_size Size of one dimension of the square grid.
// @param wall_ratio Fraction of cells to be set as walls (default: 0.1).
// @param seed Random seed for reproducibility.
// @return 2D vector of bools indicating wall positions.
std::vector<std::vector<bool>> generateRandomWallMask(int grid_size, double wall_ratio = 0.1, uint32_t seed = std::random_device{}())
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::vector<std::vector<bool>> wall_mask(grid_size, std::vector<bool>(grid_size, false));
    int wall_count = 0;

    for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
        {
            if (dist(rng) < wall_ratio)
            {
                wall_mask[i][j] = true;
                wall_count++;
            }
        }
    }

    std::cout << "Wall cells placed: " << wall_count << std::endl;
    return wall_mask;
}

// Saves a wall mask as a binary file.
// Each cell is written as a single byte (1 for wall, 0 for no wall).
// @param wall_mask The 2D boolean wall mask to save.
// @param filename The output file name.
void saveWallMaskBinary(const std::vector<std::vector<bool>> &wall_mask, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Cannot write wall mask to " << filename << std::endl;
        return;
    }

    for (const auto &row : wall_mask)
    {
        for (bool wall : row)
        {
            uint8_t value = wall ? 1 : 0;
            file.write(reinterpret_cast<char *>(&value), 1);
        }
    }
}

// Loads a wall mask from a binary file.
// @param grid_size Size of one dimension of the square grid.
// @param filename The file to load from.
// @return 2D boolean wall mask.
std::vector<std::vector<bool>> loadWallMaskBinary(int grid_size, const std::string &filename)
{
    std::vector<std::vector<bool>> wall_mask(grid_size, std::vector<bool>(grid_size, false));
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

// Applies HPP collision rules to a cell, unless it is a wall.
// If the cell is a wall, returns the cell unchanged.
// @param current_cell The cell value to process.
// @param is_wall True if the cell is a wall.
// @return The new cell value after collision.
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

// Propagates particles from the center cell to its neighbors according to HPP rules.
// Moves N, E, S, W particles from center to up, right, down, left, respectively.
// The moved particles are cleared from the center.
// @param center The center cell (by reference, particles are removed).
// @param up The cell above (by reference, receives N particles).
// @param down The cell below (by reference, receives S particles).
// @param left The cell to the left (by reference, receives W particles).
// @param right The cell to the right (by reference, receives E particles).
// Propagates particles from the center cell to its neighbors according to HPP rules.
// Moves N, E, S, W particles from center to up, right, down, left, respectively.
// The moved particles are cleared from the center.
// @param center The center cell (by reference, particles are removed).
// @param up The cell above (by reference, receives N particles).
// @param down The cell below (by reference, receives S particles).
// @param left The cell to the left (by reference, receives W particles).
// @param right The cell to the right (by reference, receives E particles).
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

// Reflects particles in a cell if it is a wall.
// Swaps N<->S and E<->W particles for wall cells; otherwise returns cell unchanged.
// @param current_cell The cell value to reflect.
// @param is_wall True if the cell is a wall.
// @return The new cell value after reflection.
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

// Inverse of the reflection operation (identical for this model).
// @param current_cell The cell value to reflect.
// @param is_wall True if the cell is a wall.
// @return The new cell value after inverse reflection.
uint8_t inverse_reflection(uint8_t current_cell, bool is_wall)
{
    return reflection(current_cell, is_wall);
}

// Inverse of the propagate operation.
// Reconstructs the center cell's particles by collecting them from the neighbors.
// @param center The center cell (by reference, particles are added).
// @param up The cell above (by reference, provides S particles).
// @param down The cell below (by reference, provides N particles).
// @param left The cell to the left (by reference, provides E particles).
// @param right The cell to the right (by reference, provides W particles).
// Inverse of the propagate operation.
// Reconstructs the center cell's particles by collecting them from the neighbors.
// @param center The center cell (by reference, particles are added).
// @param up The cell above (by reference, provides S particles).
// @param down The cell below (by reference, provides N particles).
// @param left The cell to the left (by reference, provides E particles).
// @param right The cell to the right (by reference, provides W particles).
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

// Inverse of the collision operation (identical for this model).
// @param current_cell The cell value to process.
// @param is_wall True if the cell is a wall.
// @return The new cell value after inverse collision.
uint8_t inverse_collision(uint8_t current_cell, bool is_wall)
{
    return collision(current_cell, is_wall);
}

// Main entry point for the HPP encryption and decryption simulation.
// Loads input, performs encryption and decryption using HPP CA rules, and saves results.
// @return 0 on successful completion.
int main()
{
    // ----------------------------------------------
    // Simulation of the DECRYPTION operations
    // ----------------------------------------------
    {
        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Inverse algorithm " << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Load encrypted text-file" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        std::string filename = "encrypted_message.txt";
        std::vector<uint8_t> fileData = readFileBytes(filename); // Read file as binary data

        std::cout << "Read " << fileData.size() << " bytes from " << filename << std::endl;

        // Reshape the data into a grid
        std::cout << "Reshaping data into a grid..." << std::endl;
        int grid_size;
        std::vector<std::vector<uint8_t>> grid = reshapeToMatrix(fileData, grid_size);

        // REMOVE
        std::cout << "Reshaped Matrix (" << grid_size << "x" << grid_size << "):" << std::endl;
        for (const auto &row : grid)
        {
            for (const auto &cell : row)
            {
                printBits(cell);
                std::cout << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Load wall mask... " << num_itterations << std::endl;
        std::vector<std::vector<bool>> wall_mask = loadWallMaskBinary(grid_size, "wall_mask.key");
        std::cout << "Loaded wall mask from 'wall_mask.key'" << std::endl;

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "HPP-Algorithm" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        std::cout << "Number of iterations: " << num_itterations << std::endl;
        std::cout << "Grid size:" << grid_size << "x" << grid_size << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        // Print the initial grid
        std::cout << "Initial Grid:" << std::endl;
        printGrid(grid);

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Decrypting message..." << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        for (int r = 0; r < num_itterations; ++r)
        {

            std::cout << "----------------------------------------------" << std::endl;
            std::cout << "Inverse Iteration: " << r << std::endl;
            std::cout << std::endl;
            std::cout << "Grid by iteration " << r << std::endl;
            printGrid(grid);
            std::cout << "----------------------------------------------" << std::endl;
            // ----------------------------------------------
            // Reflect the particles in the grid
            // ----------------------------------------------
            std::cout << "Reflecting particles..." << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = inverse_reflection(grid[i][j], wall_mask[i][j]);
                }
            }
            // Print the grid after reflection
            std::cout << "Grid after reflection:" << std::endl;
            printGrid(grid);

            // ----------------------------------------------
            // propagate the particles
            // ----------------------------------------------
            std::vector<std::vector<uint8_t>> propagation_grid(grid_size, std::vector<uint8_t>(grid_size, 0)); // Initialize a new grid to store the propagated values
            std::cout << "Propagating particles..." << std::endl;
            std::cout << std::endl;

            std::vector<std::vector<uint8_t>> original_grid = grid;

            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    uint8_t center_original = grid[i][j];
                    uint8_t &center = propagation_grid[i][j];
                    center = center_original & 0b11110000; // Preserve upper bits (e.g., wall)

                    uint8_t &down = original_grid[(i - 1 + grid_size) % grid_size][j];  // north particle comes from below
                    uint8_t &up = original_grid[(i + 1) % grid_size][j];                // south from above
                    uint8_t &right = original_grid[i][(j - 1 + grid_size) % grid_size]; // east from left
                    uint8_t &left = original_grid[i][(j + 1) % grid_size];              // west from right

                    inverse_propagate(center, up, down, left, right);
                }
            }

            // Copy the propagated values back to the original grid
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = propagation_grid[i][j];
                }
            }
            // Print the grid after propagation
            std::cout << "Grid after propagation:" << std::endl;
            printGrid(grid);
            std::cout << std::endl;

            // ----------------------------------------------
            // Simulate collision
            // Collision is applied to each cell in the grid
            std::cout << "Simulating collision..." << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = inverse_collision(grid[i][j], wall_mask[i][j]);
                }
            }
            // Print the grid after collision
            std::cout << "Grid after collision:" << std::endl;
            printGrid(grid);
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "End of inverse iteration " << r << std::endl;
            std::cout << "----------------------------------------------" << std::endl;
        }

        // ----------------------------------------------
        // Save final grid as ASCII text
        // ----------------------------------------------
        std::cout << "Saving final grid to ASCII text file..." << std::endl;
        std::vector<uint8_t> final_data = flattenMatrix(grid); // Flatten the 2D grid
        saveAsAsciiText(final_data, "decrypted_message.txt");  // Save as ASCII text
        std::cout << "Saved reconstructed ASCII file as 'decrypted_message.txt'" << std::endl;
    }

    return 0;
}