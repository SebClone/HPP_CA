#ifndef HPP_LIB_HPP
#define HPP_LIB_HPP

#include <vector>
#include <string>
#include <cstdint>

// Grid and file utilities
std::vector<uint8_t> readFileBytes(const std::string &filename);
std::vector<std::vector<uint8_t>> reshapeToMatrix(const std::vector<uint8_t> &data, int &grid_size);
std::vector<uint8_t> flattenMatrix(const std::vector<std::vector<uint8_t>> &matrix);
void saveAsAsciiText(const std::vector<uint8_t> &data, const std::string &filename);

// Wall mask functions
void generateRandomWallMask(std::vector<std::vector<uint8_t>> &grid, double wall_ratio = 0.1, uint32_t seed = std::random_device{}());
void saveWallMaskBinary(const std::vector<std::vector<uint8_t>> &grid, const std::string &filename);
void loadWallMaskBinary(std::vector<std::vector<uint8_t>> &grid, const std::string &filename);

// HPP Simulation functions
uint8_t collision(uint8_t current_cell);
void propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right);
uint8_t reflection(uint8_t current_cell);
uint8_t inverse_reflection(uint8_t current_cell);
void inverse_propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right);
uint8_t inverse_collision(uint8_t current_cell);

// Debug/Visualization
void printBits(uint8_t value);
void printGrid(const std::vector<std::vector<uint8_t>> &grid);

#endif // HPP_LIB_HPP