#include <iostream>
#include <cstdint> // für uint8_t
#include <bitset>  // für std::bitset

void printBits(uint8_t value)
{
    std::bitset<8> bits(value);
    std::cout << bits;
}

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
    // Grid layout (3x3): each entry is a cell with bits: xxxknosw
    uint8_t grid[3][3] = {
        {0b00000010, 0b00000000, 0b00000000},
        {0b00001000, 0b00001111, 0b00000000},
        {0b00000000, 0b00010000, 0b00000000}};

    std::cout << "Inital state of the center cell: ";
    printBits(grid[1][1]);
    std::cout << std::endl;

    // Simulate collision
    grid[1][1] = collision(grid[1][1]);
    propagate(grid[1][1], grid[0][1], grid[2][1], grid[1][0], grid[1][2]);
    uint8_t down_reflection = reflection(grid[2][1]);

    std::cout << "State of the center-middle cell after propagation: ";
    printBits(grid[1][1]);
    std::cout << std::endl;
    std::cout << "State of the up-middle cell after propagation: ";
    printBits(grid[0][1]);
    std::cout << std::endl;
    std::cout << "State of the down-middel cell after propagation: ";
    printBits(grid[2][1]);
    std::cout << std::endl;
    std::cout << "State of the left-middle cell after propagation: ";
    printBits(grid[1][0]);
    std::cout << std::endl;
    std::cout << "State of the right-middle cell after propagation: ";
    printBits(grid[1][2]);
    std::cout << std::endl;
    std::cout << "State of the left-middle cell after propagation: ";
    printBits(grid[1][2]);
    std::cout << std::endl;

    std::cout << "State of the down-middle cell after reflection: ";
    printBits(down_reflection);
    std::cout << std::endl;

    return 0;
}