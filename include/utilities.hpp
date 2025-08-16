#pragma once

#include <mpi.h>
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <random>

using Matrix = std::vector<std::vector<uint8_t>>;
using Mask = std::vector<std::vector<uint8_t>>;

// Darstellen der Ergebnisse
void printBits(uint8_t value);
void printGrid(const Matrix &grid);
void save_frame_bin(const std::vector<uint8_t> &frame, int iter);

// Binäres Einlesen
std::vector<uint8_t> readFileBytes(const std::string &filename);

// Grid-Transformationen
std::vector<std::vector<uint8_t>> reshapeToMatrix(const std::vector<uint8_t> &data, size_t &grid_size);
std::vector<uint8_t> flattenMatrix(const std::vector<std::vector<uint8_t>> &matrix);

// ASCII-Ausgabe
void saveAsAsciiText(const std::vector<uint8_t> &data, const std::string &filename);

// Wall-Masken
Mask generateRandomWallMask(int grid_size, double wall_ratio = 0.1, uint32_t seed = 0);
void saveWallMaskBinary(const Mask &wall_mask, const std::string &filename);
Mask loadWallMaskBinary(int grid_size, const std::string &filename);

// MPI-Kommunikation
void broadcastMask(Mask &mask, MPI_Comm comm);
void copy_n_bytes(const uint8_t *src, std::size_t count, uint8_t *dst);

// Gekapselte Hauptlogik der HPP-Operationen
void saveBinary(const std::vector<uint8_t> &data, const char *filename);
void saveEncryptedMeta(uint64_t originalSize, uint32_t gridSize, uint64_t &startOffset, const char *metaFilename);
bool loadEncryptedMeta(uint64_t &originalSize, uint32_t &gridSize, uint64_t &startOffset, const char *metaFilename);
