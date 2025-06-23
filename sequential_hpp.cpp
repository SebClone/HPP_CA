#include <iostream>
#include <cstdint> // für uint8_t
#include <bitset>  // für std::bitset

// ----------------------------------------------
// Custom functions for printing bits and grid
// ----------------------------------------------
void printBits(uint8_t value)
{
    std::bitset<8> bits(value);
    std::cout << bits;
}

void printGrid(uint8_t grid[3][3])
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            printBits(grid[i][j]);
            std::cout << " ";
        }
        std::cout << std::endl;
    }
}
// ----------------------------------------------

uint8_t collision(uint8_t current_cell)
{
    // Skip collision if the cell contains a wall
    if (current_cell & 0b00010000)
    {
        return current_cell;
    }

    if (!(current_cell & 0b00010000)) // No wall
    {
        if (current_cell == 0b00000101) // Cell contains a east west partical
        {
            current_cell = 0b00001010; // Cell becomes a north south partical
        }
        else if (current_cell == 0b00001010) // Cell contains a north south partical
        {
            current_cell = 0b00000101; // cell becomes a east west partical
        }
    }

    return current_cell; // Rückgabe des aktualisierten Wertes
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
    // Checks if the current cell contains a wall
    if (!(current_cell & 0b00010000))
    {
        return current_cell;
    }

    uint8_t reflected_cell = current_cell; // Initialize the reflected cell

    // Reflect the current cell by flipping the corresponding direction bits
    if (current_cell & 0b00010000) // wall bit present
    {
        if (current_cell & 0b00001000) // N
        {
            reflected_cell &= ~0b00001000;
            reflected_cell |= 0b00000010; // N → S
        }
        if (current_cell & 0b00000010) // S
        {
            reflected_cell &= ~0b00000010;
            reflected_cell |= 0b00001000; // S → N
        }
        if (current_cell & 0b00000100) // E
        {
            reflected_cell &= ~0b00000100;
            reflected_cell |= 0b00000001; // E → W
        }
        if (current_cell & 0b00000001) // W
        {
            reflected_cell &= ~0b00000001;
            reflected_cell |= 0b00000100; // W → E
        }
    }

    return reflected_cell;
}

int main()
{
    // ----------------------------------------------
    // Initialization of the grid
    // Grid layout (3x3): each entry is a cell with bits: xxxknesw
    // ----------------------------------------------
    std::cout << "Particle Bit" << std::endl;
    std::cout << "k: wall, n: north, e: east, s: south, w: west" << std::endl;
    std::cout << "Example: xxxkNESW" << std::endl;

    uint8_t grid[3][3] = {
        {0b00000010, 0b00000000, 0b00001000},
        {0b00001000, 0b00001111, 0b00000000},
        {0b00000101, 0b00010000, 0b00001010}};

    // Print the initial grid
    std::cout << "Initial Grid:" << std::endl;
    printGrid(grid);

    // ----------------------------------------------
    // Simulate collision
    // Collision is applied to each cell in the grid

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            grid[i][j] = collision(grid[i][j]);
        }
    }
    // Print the grid after collision
    std::cout << "Grid after collision:" << std::endl;
    printGrid(grid);

    // ----------------------------------------------
    // propagate the particles
    // ----------------------------------------------
    uint8_t propagation_grid[3][3] = {0}; // Initialize a new grid to store the propagated values

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // Check wehter there is a next cell to propagate to
            uint8_t &center = grid[i][j]; // current center cell

            // Toroidal neighbor assignments
            // Using + 3 and % 3 to move in a toroidal manner
            uint8_t &up = propagation_grid[(i - 1 + 3) % 3][j];
            uint8_t &down = propagation_grid[(i + 1) % 3][j];
            uint8_t &left = propagation_grid[i][(j - 1 + 3) % 3];
            uint8_t &right = propagation_grid[i][(j + 1) % 3];

            uint8_t temp = center;
            propagate(temp, up, down, left, right);
            propagation_grid[i][j] |= temp; // Not prpagatet particals will remain in the gird and not be deleted
        }
    }

    // Copy the propagated values back to the original grid
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            grid[i][j] = propagation_grid[i][j];
        }
    }
    // Print the grid after propagation
    std::cout << "Grid after propagation:" << std::endl;
    printGrid(grid);

    // ----------------------------------------------
    // Reflect the particles in the grid
    // ----------------------------------------------
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            grid[i][j] = reflection(grid[i][j]);
        }
    }
    // Print the grid after reflection
    std::cout << "Grid after reflection:" << std::endl;
    printGrid(grid);

    return 0;
}