#include <iostream>
#include <cstdint> // für uint8_t
#include <bitset>  // für std::bitset
#include <cstring> // for memcpy
#include <fstream>
#include <vector>

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
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "Load Text-File" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::string filename = "message.txt";
    std::vector<uint8_t> fileData = readFileBytes(filename); // Read file as binary data

    std::cout << "Read " << fileData.size() << " bytes from " << filename << std::endl;

    // Print first few bytes
    for (size_t i = 0; i < std::min<size_t>(10, fileData.size()); ++i)
    {
        std::cout << std::hex << static_cast<int>(fileData[i]) << " ";
    }
    std::cout << std::endl;

    // Reshape the data into a grid
    std::cout << "Reshaping data into a grid..." << std::endl;
    int grid_size;
    std::vector<std::vector<uint8_t>> grid = reshapeToMatrix(fileData, grid_size);

    std::cout << "Reshaped Matrix (" << grid_size << "x" << grid_size << "):" << std::endl;
    for (const auto &row : grid) // Remove
    {
        for (const auto &cell : row)
        {
            printBits(cell);
            std::cout << " ";
        }
        std::cout << std::endl;
    }

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
    // Simulation of the inverse operations
    // ----------------------------------------------
    for (int r = 0; r < num_itterations; ++r)
    {
        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Inverse algorithm " << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

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
    saveAsAsciiText(final_data, "output.txt");             // Save as ASCII text
    std::cout << "Saved reconstructed ASCII file as 'output.txt'" << std::endl;

    return 0;
}