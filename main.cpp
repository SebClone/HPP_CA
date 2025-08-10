#include "hpp_encryptor.hpp"

constexpr int NUM_ITERATIONS = 1000;

int main()
{
    bool doEncrypt       = true;   
    bool doDecrypt       = true; 
    bool visualize       = false;  
    bool printEncrypted  = true;
    bool printDecrypted  = true;

    if (!doEncrypt && !doDecrypt) 
    {
        std::cerr << "Please enable encrypt or decrypt in the flags.\n";
        return 1;
    }

    // ----------------------------------------------
    // ENCRYPTION 
    // ----------------------------------------------
    if (doEncrypt) 
    {
        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Load Text-File" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        std::string filename = "message.txt";
        std::vector<uint8_t> fileData = readFileBytes(filename); // Read file as binary data

        std::cout << "Read " << fileData.size() << " bytes from " << filename << std::endl;

        std::cout << "Reshaping data into a grid..." << std::endl;
        size_t grid_size;
        std::vector<std::vector<uint8_t>> grid = reshapeToMatrix(fileData, grid_size);

        if (visualize)
        {
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
        }

        std::cout << "Generating random wall mask..." << std::endl;
        std::cout << "This will set ~10% of the grid cells to walls." << std::endl;
        Mask wall_mask = generateRandomWallMask(grid_size, 0.1);
        std::cout << "Saving wall mask to binary file..." << std::endl;
        saveWallMaskBinary(wall_mask, "wall_mask.key");
        std::cout << "Saved wall mask to 'wall_mask.key'" << std::endl;

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "HPP-Algorithm" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        std::cout << "Number of iterations: " << NUM_ITERATIONS << std::endl;
        std::cout << "Grid size:" << grid_size << "x" << grid_size << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        if (visualize)
        {
            std::cout << "Initial Grid:" << std::endl;
            printGrid(grid);
        }

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Encrypting message..." << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        for (int r = 0; r < NUM_ITERATIONS; ++r)
        {
            if (visualize)
            {
                std::cout << "----------------------------------------------" << std::endl;
                std::cout << "Iteration: " << r << std::endl;
                std::cout << std::endl;
                std::cout << "Grid by iteration " << r << std::endl;
                printGrid(grid);
                std::cout << "----------------------------------------------" << std::endl;
            
                std::cout << "Simulating collision..." << std::endl;
                std::cout << std::endl;
            }

            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = collision(grid[i][j], wall_mask[i][j]);
                }
            }
            
            if (visualize)
            {
                std::cout << "Grid after collision:" << std::endl;
                printGrid(grid);
                std::cout << std::endl;
            }

            std::vector<std::vector<uint8_t>> propagation_grid(grid_size, std::vector<uint8_t>(grid_size, 0)); // Initialize a new grid to store the propagated values
            
            if (visualize)
            {
                std::cout << "Propagating particles..." << std::endl;
                std::cout << std::endl;
            }
    
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    // Check whether there is a next cell to propagate to
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

            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = propagation_grid[i][j];
                }
            }
            
            if (visualize)
            {
                std::cout << "Grid after propagation:" << std::endl;
                printGrid(grid);
                std::cout << std::endl;
            }

            if (visualize)
            {
                std::cout << "Reflecting particles..." << std::endl;
                std::cout << std::endl;
            }
        
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = reflection(grid[i][j], wall_mask[i][j]);
                }
            }
            
            if (visualize)
            {
                std::cout << "Grid after reflection:" << std::endl;
                printGrid(grid);
                std::cout << std::endl;
                std::cout << "End of iteration " << r << std::endl;
                std::cout << "----------------------------------------------" << std::endl;
            }
        }

        std::cout << "Saving final grid to ASCII text file..." << std::endl;
        std::vector<uint8_t> encrypted_data = flattenMatrix(grid);

        saveAsAsciiText(encrypted_data, "encrypted_message.txt"); 
        std::cout << "Saved reconstructed ASCII file as 'encrypted_message.txt'" << std::endl;

        if (printEncrypted)
        {
            std::cout << "Encrypted message (as ASCII):" << std::endl;
            for (auto c : encrypted_data) {
                std::cout << static_cast<char>(c);
            }
            std::cout << std::endl;
        }
    }


    // ----------------------------------------------
    // DECRYPTION
    // ----------------------------------------------
    if (doDecrypt) 
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

        std::cout << "Reshaping data into a grid..." << std::endl;
        size_t grid_size;
        std::vector<std::vector<uint8_t>> grid = reshapeToMatrix(fileData, grid_size);

        if (visualize)
        {
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
        }

        std::cout << "Load wall mask... " << NUM_ITERATIONS << std::endl;
        Mask wall_mask = loadWallMaskBinary(grid_size, "wall_mask.key");
        std::cout << "Loaded wall mask from 'wall_mask.key'" << std::endl;

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "HPP-Algorithm" << std::endl;
        std::cout << "----------------------------------------------" << std::endl;

        std::cout << "Number of iterations: " << NUM_ITERATIONS << std::endl;
        std::cout << "Grid size:" << grid_size << "x" << grid_size << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        
        if (visualize)
        {
            std::cout << "Initial Grid:" << std::endl;
            printGrid(grid);
        }

        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Decrypting message..." << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        for (int r = 0; r < NUM_ITERATIONS; ++r)
        {
            if (visualize)
            {
                std::cout << "----------------------------------------------" << std::endl;
                std::cout << "Inverse Iteration: " << r << std::endl;
                std::cout << std::endl;
                std::cout << "Grid by iteration " << r << std::endl;
                printGrid(grid);
                std::cout << "----------------------------------------------" << std::endl;
            }
            
            if (visualize)
            {
                std::cout << "Reflecting particles..." << std::endl;
                std::cout << std::endl;
            }

            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = inverse_reflection(grid[i][j], wall_mask[i][j]);
                }
            }
            
            if (visualize)
            {
                std::cout << "Grid after reflection:" << std::endl;
                printGrid(grid);
            }

            std::vector<std::vector<uint8_t>> propagation_grid(grid_size, std::vector<uint8_t>(grid_size, 0)); // Initialize a new grid to store the propagated values

            if (visualize)
            {
                std::cout << "Propagating particles..." << std::endl;
                std::cout << std::endl;
            }

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

            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = propagation_grid[i][j];
                }
            }
            
            if (visualize)
            {
                std::cout << "Grid after propagation:" << std::endl;
                printGrid(grid);
                std::cout << std::endl;
            }

            if (visualize)
            {
                std::cout << "Simulating collision..." << std::endl;
                std::cout << std::endl;
            }
            
            for (int i = 0; i < grid_size; ++i)
            {
                for (int j = 0; j < grid_size; ++j)
                {
                    grid[i][j] = inverse_collision(grid[i][j], wall_mask[i][j]);
                }
            }
            
            if (visualize)
            {
                std::cout << "Grid after collision:" << std::endl;
                printGrid(grid);
                std::cout << std::endl;
                std::cout << std::endl;
                std::cout << "End of inverse iteration " << r << std::endl;
                std::cout << "----------------------------------------------" << std::endl;
            }
        }

        std::cout << "Saving final grid to ASCII text file..." << std::endl;
        std::vector<uint8_t> final_data = flattenMatrix(grid); // Flatten the 2D grid
        saveAsAsciiText(final_data, "decrypted_message.txt");  // Save as ASCII text
        std::cout << "Saved reconstructed ASCII file as 'decrypted_message.txt'" << std::endl;

        if (printDecrypted)
        {
            std::cout << "Decrypted message (as ASCII):" << std::endl;
            for (auto c : final_data) {
                std::cout << static_cast<char>(c);
            }
            std::cout << std::endl;
        }
    }

    return 0;
}