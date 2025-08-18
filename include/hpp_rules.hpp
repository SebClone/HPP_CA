#pragma once

#include <cstdint>
#include <vector>
#include <cstddef>

// Particle-Kollision & -Propagation
uint8_t collision(uint8_t current_cell, bool is_wall);
void propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right);
uint8_t reflection(uint8_t current_cell, bool is_wall);
uint8_t inverse_reflection(uint8_t current_cell, bool is_wall);
void inverse_propagate(uint8_t &center, uint8_t &up, uint8_t &down, uint8_t &left, uint8_t &right);
uint8_t inverse_collision(uint8_t current_cell, bool is_wall);

//Gesamtlogik der HPP-Operationen
using Matrix = std::vector<std::vector<uint8_t>>;
using Mask = std::vector<std::vector<uint8_t>>;

uint8_t applyRules(
    const Matrix& active_grid,
    const Mask&   wall_mask,
    bool          doEncrypt,
    int           i,
    int           j,
    int           offset_rows);


// Portables "restrict"-Makro: Compiler erhält einen Alias-Hinweis für bessere Vektorisierung.
#if !defined(FAST_RESTRICT)
  #if defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
    #define FAST_RESTRICT __restrict
  #else
    #define FAST_RESTRICT
  #endif
#endif

template<bool ENCRYPT>
inline uint8_t applyRules_fast(
    const uint8_t* FAST_RESTRICT G,
    int rowStride,                 // lokale Zeilenbreite inkl. Halos (W = local_cols+2)
    int i, int j,                  // lokale Indizes inkl. Halos (i ∈ [1..local_rows], j ∈ [1..local_cols])
    const uint8_t* FAST_RESTRICT wrow,     // BEGINN der globalen Maskenzeile
    const uint8_t* FAST_RESTRICT wrow_up,  // BEGINN der globalen Maskenzeile oben
    const uint8_t* FAST_RESTRICT wrow_dn,  // BEGINN der globalen Maskenzeile unten
    int Nmask,                     // globale Maskenbreite (= grid_size)
    int j0_global                  // globale Spalte, die lokal bei j==1 liegt (= offset_cols)
);