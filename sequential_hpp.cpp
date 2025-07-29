#include <iostream>
#include <cstdint> // für uint8_t
#include <bitset>  // für std::bitset
#include <cstring> // for memcpy
#include <fstream>
#include <vector>
#include <random>

const int num_itterations = 1; // Number of iterations for the simulation

// ----------------------------------------------
// Custom functions for printing bits and grid
// ----------------------------------------------
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

// Reads a file as an binary and returns its contents as a vector of bytes
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

/**
 * @brief Randomly generates wall cells in the grid by setting the wall bit (bit 4) in each cell
 *        with probability wall_ratio. Uses a random number generator seeded with 'seed'.
 *
 * @param grid The 2D grid whose cells will be assigned wall bits.
 * @param wall_ratio The probability (0.0 to 1.0) that any given cell will become a wall.
 * @param seed Random seed for reproducibility.
 *
 * The function iterates over every cell in the grid. For each cell, it generates a random number in [0, 1).
 * If the random number is less than wall_ratio, the cell is turned into a wall by setting bit 4 (0b00010000).
 * The function also counts and prints the number of wall cells placed.
 */
void generateRandomWallMask(std::vector<std::vector<uint8_t>> &grid, double wall_ratio = 0.1, uint32_t seed = std::random_device{}())
{
    std::mt19937 rng(seed);                                // Mersenne Twister RNG seeded with 'seed'
    std::uniform_real_distribution<double> dist(0.0, 1.0); // Uniform distribution [0,1)

    int wall_count = 0; // Counter for number of wall cells placed
    for (auto &row : grid)
    {
        for (auto &cell : row)
        {
            // For each cell, randomly decide if it becomes a wall
            if (dist(rng) < wall_ratio)
            {
                cell |= 0b00010000; // Set wall bit (bit 4)
                wall_count++;
            }
        }
    }

    std::cout << "Wall cells placed: " << wall_count << std::endl;
}

/**
 * @brief Saves the wall mask of the grid to a binary file.
 *
 * @param grid The 2D grid whose wall mask will be saved.
 * @param filename The name of the file to save to.
 *
 * The function iterates over every cell in the grid, extracting only the wall bit (bit 4).
 * For each cell, it writes 1 (if the wall bit is set) or 0 (if not) as a single byte to the file.
 * The result is a binary file containing only wall information for each cell, in row-major order.
 */
void saveWallMaskBinary(const std::vector<std::vector<uint8_t>> &grid, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Cannot write wall mask to " << filename << std::endl;
        return;
    }

    for (const auto &row : grid)
    {
        for (uint8_t cell : row)
        {
            // Extract wall bit: write 1 if wall, 0 otherwise
            uint8_t wall = (cell & 0b00010000) ? 1 : 0;
            file.write(reinterpret_cast<char *>(&wall), 1);
        }
    }
}

/**
 * @brief Loads the wall mask from a binary file and applies it to the grid.
 *
 * @param grid The 2D grid to which the wall mask will be applied.
 * @param filename The name of the binary file to load from.
 *
 * The function reads one byte per cell from the file, in row-major order.
 * If the byte is 1, the wall bit (bit 4) is set in the corresponding cell.
 * If the byte is 0, the wall bit is cleared.
 * This overwrites the wall bit in each cell but leaves other bits unchanged.
 */
void loadWallMaskBinary(std::vector<std::vector<uint8_t>> &grid, const std::string &filename)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Cannot read wall mask from " << filename << std::endl;
        return;
    }

    for (auto &row : grid)
    {
        for (auto &cell : row)
        {
            char wall;
            file.read(&wall, 1); // Read one byte for wall info
            if (wall == 1)
                cell |= 0b00010000; // Set wall bit (bit 4)
            else
                cell &= ~0b00010000; // Clear wall bit (bit 4)
        }
    }
}

// ----------------------------------------------
// Forward declarations of functions
// ----------------------------------------------
uint8_t collision(uint8_t current_cell)
{
    // Skip collision if the cell contains a wall
    if (current_cell & 0b00010000)
    {
        return current_cell;
    }

    uint8_t particles = current_cell & 0b00001111; // Mask out only NESW bits

    if (particles == 0b00000101) // E + W
    {
        current_cell &= 0b11110000; // Clear NESW
        current_cell |= 0b00001010; // Set N + S
    }
    else if (particles == 0b00001010) // N + S
    {
        current_cell &= 0b11110000; // Clear NESW
        current_cell |= 0b00000101; // Set E + W
    }

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
        right |= 0b00000100;   // up cell becomes a east partical
        center &= ~0b00000100; // clear the east partical from center
    }
    // Propagation if its an south partical
    if (center & 0b00000010)
    {
        down |= 0b00000010;    // up cell becomes a south partical
        center &= ~0b00000010; // clear the south partical from center
    }
    // Propagation if its an south partical
    if (center & 0b00000001)
    {
        left |= 0b00000001;    // up cell becomes a west partical
        center &= ~0b00000001; // clear the west partical from center
    }
}

uint8_t reflection(uint8_t current_cell)
{
    // Return early if there's no wall
    if (!(current_cell & 0b00010000))
    {
        return current_cell;
    }

    // Extract NESW bits only
    uint8_t particles = current_cell & 0b00001111;
    uint8_t reflected_particles = 0;

    // Reflect each particle direction
    if (particles & 0b00001000) // N → S
        reflected_particles |= 0b00000010;

    if (particles & 0b00000010) // S → N
        reflected_particles |= 0b00001000;

    if (particles & 0b00000100) // E → W
        reflected_particles |= 0b00000001;

    if (particles & 0b00000001) // W → E
        reflected_particles |= 0b00000100;

    // Combine reflected particles with the original non-NESW bits (wall, etc.)
    uint8_t non_particle_bits = current_cell & 0b11110000;
    return non_particle_bits | reflected_particles;
}

// ----------------------------------------------
// Backward declarations of functions
// ----------------------------------------------

uint8_t inverse_reflection(uint8_t current_cell)
{
    return reflection(current_cell); // Inverse reflection is the same as reflection
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

uint8_t inverse_collision(uint8_t current_cell)
{
    return collision(current_cell); // Inverse collision is the same as collision
}

int main()
{
    // ----------------------------------------------
    // Simulation of the ENCRYPTION operations
    // ----------------------------------------------
    {
        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Load Text-File" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        std::string filename = "message.txt";
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

        // Generate a random wall mask and save it to a binary file
        std::cout << "Generating random wall mask..." << std::endl;
        std::cout << "This will set ~10% of the grid cells to walls." << std::endl;
        generateRandomWallMask(grid, 0.1); // ~10% of grid cells become walls
        std::cout << "Saving wall mask to binary file..." << std::endl;
        saveWallMaskBinary(grid, "wall_mask.key");
        std::cout << "Saved wall mask to 'wall_mask.key'" << std::endl;

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "HPP-Algorithm" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        std::cout << "Number of iterations: " << num_itterations << std::endl;
        std::cout << "Grid size:" << grid_size << "x" << grid_size << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        // ----------------------------------------------
        // Initialization of the grid
        // Grid layout (grid_sizexgrid_size): each entry is a cell with bits: xxxknesw
        // ----------------------------------------------

        // Print the initial grid
        std::cout << "Initial Grid:" << std::endl;
        printGrid(grid);

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Encrypting message..." << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        for (int r = 0; r < num_itterations; ++r)
        {
            std::cout << "----------------------------------------------" << std::endl;
            std::cout << "Iteration: " << r << std::endl;
            std::cout << std::endl;
            std::cout << "Grid by iteration " << r << std::endl;
            printGrid(grid);
            std::cout << "----------------------------------------------" << std::endl;
            // ----------------------------------------------
            // Simulate collision
            // Collision is applied to each cell in the grid
            std::cout << "Simulating collision..." << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = collision(grid[i][j]);
                }
            }
            // Print the grid after collision
            std::cout << "Grid after collision:" << std::endl;
            printGrid(grid);
            std::cout << std::endl;

            // ----------------------------------------------
            // propagate the particles
            // ----------------------------------------------
            std::vector<std::vector<uint8_t>> propagation_grid(grid_size, std::vector<uint8_t>(grid_size, 0)); // Initialize a new grid to store the propagated values
            std::cout << "Propagating particles..." << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    // Check wehter there is a next cell to propagate to
                    uint8_t &center = grid[i][j]; // current center cell

                    // Toroidal neighbor assignments
                    // Using + grid_size and % grid_size to move in a toroidal manner
                    uint8_t &up = propagation_grid[(i - 1 + grid_size) % grid_size][j];
                    uint8_t &down = propagation_grid[(i + 1) % grid_size][j];
                    uint8_t &left = propagation_grid[i][(j - 1 + grid_size) % grid_size];
                    uint8_t &right = propagation_grid[i][(j + 1) % grid_size];

                    uint8_t temp = center;
                    propagate(temp, up, down, left, right);
                    propagation_grid[i][j] |= temp; // Not prpagatet particals will remain in the gird and not be deleted
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
            // Reflect the particles in the grid
            // ----------------------------------------------
            std::cout << "Reflecting particles..." << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = reflection(grid[i][j]);
                }
            }
            // Print the grid after reflection
            std::cout << "Grid after reflection:" << std::endl;
            printGrid(grid);
            std::cout << std::endl;
            std::cout << "End of iteration " << r << std::endl;
            std::cout << "----------------------------------------------" << std::endl;
        }

        // ----------------------------------------------
        // Save encrypted grid as ASCII text
        // ----------------------------------------------
        std::cout << "Saving final grid to ASCII text file..." << std::endl;
        std::vector<uint8_t> encrypted_data = flattenMatrix(grid); // Flatten the 2D grid

        // Strip wall bit (bit 4) before saving
        for (uint8_t &byte : encrypted_data)
        {
            byte &= 0b11101111; // Clear bit 4 (wall bit)
        }

        saveAsAsciiText(encrypted_data, "encrypted_message.txt"); // Save as ASCII text
        std::cout << "Saved reconstructed ASCII file as 'encrypted_message.txt'" << std::endl;
    }

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
        loadWallMaskBinary(grid, "wall_mask.key");
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
                    grid[i][j] = inverse_reflection(grid[i][j]);
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
                    grid[i][j] = inverse_collision(grid[i][j]);
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