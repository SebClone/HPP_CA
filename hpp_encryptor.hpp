// Enthält die Deklarationen aller Funktionen und alle include statements

#ifndef HPP_ENCRYPTOR_HPP
#define HPP_ENCRYPTOR_HPP

#include <iostream>
#include <cstdint> // für uint8_t
#include <bitset>  // für std::bitset
#include <cstring> // for memcpy
#include <fstream>
#include <vector>
#include <random>
#include <mpi.h>


//Darstellen der Ergebnisse
void printBits(uint8_t value);
void printGrid(const std::vector<std::vector<uint8_t>> &grid);

// Binäres Einlesen
std::vector<uint8_t> readFileBytes(const std::string &filename);

// Grid-Transformationen
std::vector<std::vector<uint8_t>> reshapeToMatrix(const std::vector<uint8_t> &data, size_t &grid_size);
std::vector<uint8_t> flattenMatrix(const std::vector<std::vector<uint8_t>> &matrix);

// ASCII-Ausgabe
void saveAsAsciiText(const std::vector<uint8_t> &data, const std::string &filename);

// Wall-Masken
using Mask = std::vector<std::vector<uint8_t>>;
Mask generateRandomWallMask(int grid_size, double wall_ratio = 0.1, uint32_t seed = std::random_device{}());
void saveWallMaskBinary(const Mask &wall_mask, const std::string &filename);
Mask loadWallMaskBinary(int grid_size, const std::string &filename);

// Particle-Kollision & -Propagation
uint8_t collision(uint8_t current_cell, bool is_wall);
void propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right);
uint8_t reflection(uint8_t current_cell, bool is_wall);
uint8_t inverse_reflection(uint8_t current_cell, bool is_wall);
void inverse_propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right);
uint8_t inverse_collision(uint8_t current_cell, bool is_wall);

// MPI-Kommunikation
using Matrix = std::vector<std::vector<uint8_t>>;
void broadcastMask(Mask& mask, MPI_Comm comm);
void copy_n_bytes(const uint8_t* src, std::size_t count, uint8_t* dst);
uint8_t applyRules(Matrix& grid, const Mask& wall_mask,bool forward, int row,int col);

#endif 

