#pragma once

#include <mpi.h>
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

// Beschreibt den lokalen Block dieses Ranks
struct RowDist
{
    int grid_size = 0;   // N (Spalten)
    int local_rows = 0;  // #Zeilen dieses Ranks (ohne Halos)
    int offset_rows = 0; // Startzeile im globalen Raster
};

// --- Plain/Cipher ---

// Liest Klartext in den lokalen Block (Padding: restliche Bytes = 0)
void parallel_read_plain_chunk(const std::string &path,
                               const RowDist &dist,
                               std::uint64_t original_size,
                               std::uint64_t start_offset,
                               std::vector<std::uint8_t> &out, // size == local_rows*N
                               bool atomic_io,
                               MPI_Comm comm);

// Liest Cipher (verschlüsseltes Raster, exakt N*N Bytes)
void parallel_read_cipher_chunk(const std::string &path,
                                const RowDist &dist,
                                std::size_t padded_size_bytes, // N*N
                                std::vector<std::uint8_t> &out,
                                bool atomic_io,
                                MPI_Comm comm);

// Schreibt Cipher (N*N) – Rank 0 truncatet/expandiert zuerst.
void parallel_write_cipher_chunk(const std::string &path,
                                 const RowDist &dist,
                                 const std::vector<std::uint8_t> &data,
                                 std::size_t padded_size_bytes,
                                 bool atomic_io,
                                 MPI_Comm comm);

// Schreibt Klartext (auf originalSize gekürzt)
void parallel_write_plain_trimmed(const std::string &path,
                                  const RowDist &dist,
                                  const std::vector<std::uint8_t> &data,
                                  std::uint64_t original_size,
                                  std::uint64_t start_offset,
                                  bool atomic_io,
                                  MPI_Comm comm);

// Optional: Frame-Dump (N*N, für Debug/Visualisierung)
void dump_frame_parallel(const std::string &path,
                         const RowDist &dist,
                         const std::vector<std::uint8_t> &data,
                         std::size_t padded_size_bytes,
                         bool atomic_io,
                         MPI_Comm comm);

// --- Meta (klein; Rank 0 I/O) ---
bool read_meta_rank0(const std::string &meta_path,
                     std::uint64_t &original_size,
                     std::uint32_t &grid_size_n,
                     std::uint64_t &start_offset);

bool write_meta_rank0(const std::string &meta_path,
                      std::uint64_t original_size,
                      std::uint32_t grid_size_n,
                      std::uint64_t start_offset);
