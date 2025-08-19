#include "hpp_rules.hpp"

/**
 * @brief Applies the HPP rules for particle collision to a given cell.
 *
 * Takes the current cell and corresponding the wall mask to determine wehter to apply the collision rules.
 * If the cell is a wall, it returns the current cell unchanged.
 * If the cell is not a wall, it checks for specific particle configurations and applies the collision rules accordingly.
 * E.g. xxxx1010 -> 00001010 (north-south particle collide -> east-west particle)
 * @param current_cell The current cell value (8 bits).
 * @param is_wall A boolean indicating whether the current cell is a wall (true) or not (false).
 * @return The updated cell value after applying the collision rules.
 */
uint8_t collision(uint8_t current_cell, bool is_wall)
{
    if (is_wall)
        return current_cell;

    uint8_t particles = current_cell & 0b00001111;

    if (particles == 0b00000101)
        return (current_cell & 0b11110000) | 0b00001010;
    else if (particles == 0b00001010)
        return (current_cell & 0b11110000) | 0b00000101;

    return current_cell;
}

/**
 * @brief Propagates particles in the HPP model. SEQUENTIAL
 *
 * Function is best used in a sequential context. For parallel processing, use the `applyRules_fast` function.
 * Takes the current cell and its neighbors, and propagates particles according to the HPP rules.
 * NOTE: For more details see the HPP_states.mdm
 * @param center The current cell value (8 bits).
 * @param up The value of the cell above the current cell.
 * @param down The value of the cell below the current cell.
 * @param left The value of the cell to the left of the current cell.
 * @param right The value of the cell to the right of the current cell.
 */
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
        right |= 0b00000100;   // right cell becomes a east partical
        center &= ~0b00000100; // clear the east partical from center
    }
    // Propagation if its an south partical
    if (center & 0b00000010)
    {
        down |= 0b00000010;    // down cell becomes a south partical
        center &= ~0b00000010; // clear the south partical from center
    }
    // Propagation if its an south partical
    if (center & 0b00000001)
    {
        left |= 0b00000001;    // left cell becomes a west partical
        center &= ~0b00000001; // clear the west partical from center
    }
}

/**
 * @brief Applies the HPP rules for particle reflection.
 *
 * Takes the current cell and the wall mask to determine whether to apply the reflection rules.
 * If the cell is a wall, it reflects the particles according to the HPP rules.
 * particle directions are rotated by 180 degrees e.g. 00001000 -> 00000010 (north particle reflect -> south particle)
 * @param current_cell The current cell value (8 bits).
 * @param is_wall A boolean indicating whether the current cell is a wall (true) or not (false).
 * @return The updated cell value after applying the reflection rules.
 */
uint8_t reflection(uint8_t current_cell, bool is_wall)
{
    if (!is_wall)
        return current_cell;

    uint8_t particles = current_cell & 0b00001111;
    uint8_t reflected = 0;
    if (particles & 0b00001000)
        reflected |= 0b00000010;
    if (particles & 0b00000010)
        reflected |= 0b00001000;
    if (particles & 0b00000100)
        reflected |= 0b00000001;
    if (particles & 0b00000001)
        reflected |= 0b00000100;

    return (current_cell & 0b11110000) | reflected;
}

/**
 * @brief Applies the inverse reflection rules for particles in the HPP model.
 *
 * Takes the current cell and the wall mask to determine whether to apply the inverse reflection rules.
 * If the cell is a wall, it reflects the particles according to the HPP rules in reverse.
 * particle directions are rotated by 180 degrees e.g. 00001000 -> 00000010 (north particle reflect -> south particle)
 * NOTE: Due to the nature of the HPP model, the inverse reflection is the same as the reflection.
 * @param current_cell The current cell value (8 bits).
 * @param is_wall A boolean indicating whether the current cell is a wall (true) or not (false).
 * @return The updated cell value after applying the inverse reflection rules.
 */
uint8_t inverse_reflection(uint8_t current_cell, bool is_wall)
{
    return reflection(current_cell, is_wall);
}

/**
 * @brief Applies the inverse propagation rules for particles in the HPP model.
 *
 * Takes the current cell and its neighbors, and propagates particles according to the HPP rules in reverse.
 * This is used in decryption to reconstruct the original particle configuration.
 * North particles are taken from the UP cell, South particles from the DOWN cell,
 * East particles from the RIGHT cell, and West particles from the LEFT cell.
 * @param center The current cell value (8 bits).
 * @param up The value of the cell above the current cell.
 * @param down The value of the cell below the current cell.
 * @param left The value of the cell to the left of the current cell.
 * @param right The value of the cell to the right of the current cell.
 */
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

/**
 * @brief Applies the inverse collision rules for particles in the HPP model.
 *
 * Takes the current cell and the wall mask to determine whether to apply the inverse collision rules.
 * If the cell is a wall, it returns the current cell unchanged.
 * If the cell is not a wall, it checks for specific particle configurations and applies the inverse collision rules accordingly.
 * E.g. 00001010 -> 00001010 (east-west particle collide -> north-south particle)
 * NOTE: Due to the nature of the HPP model, the inverse collision is the same as the collision.
 * @param current_cell The current cell value (8 bits).
 * @param is_wall A boolean indicating whether the current cell is a wall (true) or not (false).
 * @return The updated cell value after applying the inverse collision rules.
 */
uint8_t inverse_collision(uint8_t current_cell, bool is_wall)
{
    return collision(current_cell, is_wall);
}

// Gekapselte Hauptlogik der HPP-Operationen

/**
 * @brief Applies the HPP rules for particle propagation, collision, and reflection.
 *
 * Can be used for both encryption and decryption.
 * Takes the active grid, wall mask, and current cell coordinates to apply the HPP rules.
 * @param active_grid The current state of the grid (2D vector of bytes).
 * @param wall_mask The wall mask indicating where walls are present in the grid.
 * @param doEncrypt A boolean indicating whether to encrypt (true) or decrypt (false).
 * @param i The local row index in the subgrid (with halo): 1...local_rows
 * @param j The column index: 0..N-1
 * @param offset_rows The global starting row of this rank (without halo).
 * @return
 */
uint8_t applyRules(
    const Matrix &active_grid,
    const Mask &wall_mask,
    bool doEncrypt,  // true = forward (Encrypt), false = backward (Decrypt)
    int i,           // local rows in the sub-lattice (with Halo): 1..local_rows
    int j,           // cloumn: 0..N-1
    int offset_rows) // global starting row of this rank (without halo)
{
    const int N = static_cast<int>(wall_mask.size()); // squared lattice size

    // Global coordinates of the target cell (rows via rank, columns locally in the torus)
    const int gr = (offset_rows + (i - 1) + N) % N;
    const int gc = j;

    // Neigbouring indices (local and global)
    const int iu = i - 1;           // Neighbour from above  (recived from Halo/Inside)
    const int id = i + 1;           // Neighbour from blow (recived from Halo/Inside)
    const int jl = (j - 1 + N) % N; // Neighbour to the left  (torus)
    const int jr = (j + 1) % N;     // Nieghbour to the right (torus)

    // Load rows (read only!)
    const uint8_t c = active_grid[i][j];
    const uint8_t up = active_grid[iu][j];
    const uint8_t dn = active_grid[id][j];
    const uint8_t lf = active_grid[i][jl];
    const uint8_t rt = active_grid[i][jr];

    // Wall-Flags (Mask is uint8_t → nonzero == true)
    const bool w_c = wall_mask[gr][gc] != 0;
    const bool w_up = wall_mask[(gr - 1 + N) % N][gc] != 0;
    const bool w_dn = wall_mask[(gr + 1) % N][gc] != 0;
    const bool w_lf = wall_mask[gr][jl] != 0;
    const bool w_rt = wall_mask[gr][jr] != 0;

    if (doEncrypt)
    {
        // ---------- Forward: collision -> propagate (collect in traget row) -> reflection ----------

        // 1) Apply collision on center + neighbours
        const uint8_t c_col = collision(c, w_c);
        const uint8_t up_col = collision(up, w_up);
        const uint8_t dn_col = collision(dn, w_dn);
        const uint8_t lf_col = collision(lf, w_lf);
        const uint8_t rt_col = collision(rt, w_rt);

        // 2) Aplly "propagation" to these cells:
        //    Keep upper 4 bits from center, collect direction bits from collided neighbors
        uint8_t next = static_cast<uint8_t>(c_col & 0b11110000);

        if (up_col & 0b00000010)
            next |= 0b00000010; // From ABOVE comes their SOUTH
        if (dn_col & 0b00001000)
            next |= 0b00001000; // From BELOW comes their NORTH
        if (lf_col & 0b00000100)
            next |= 0b00000100; // From LEFT comes their EAST
        if (rt_col & 0b00000001)
            next |= 0b00000001; // From RIGHT comes their WEST

        // 3) Apply reflection to target cell
        return reflection(next, w_c);
    }
    else
    {
        // ---------- BACKWARD/INVERSE: inverse_reflection -> inverse_propagate (collect) -> inverse_collision ----------

        // 1) Apply inverse_reflection to center cell + neighbours
        const uint8_t c_ref = inverse_reflection(c, w_c);
        const uint8_t up_ref = inverse_reflection(up, w_up);
        const uint8_t dn_ref = inverse_reflection(dn, w_dn);
        const uint8_t lf_ref = inverse_reflection(lf, w_lf);
        const uint8_t rt_ref = inverse_reflection(rt, w_rt);

        // 2) Apply inverse_propagate to cell (corresponds to grid-wide variant):
        uint8_t prev = static_cast<uint8_t>(c_ref & 0b11110000);

        if (up_ref & 0b00001000)
            prev |= 0b00001000; // NORTH comes from ABOVE
        if (dn_ref & 0b00000010)
            prev |= 0b00000010; // SOUTH comes from BELOW
        if (rt_ref & 0b00000100)
            prev |= 0b00000100; // EAST comes from RIGHT
        if (lf_ref & 0b00000001)
            prev |= 0b00000001; // WEST comes from LEFT

        // 3) Apply inverse_collision to target cell
        return inverse_collision(prev, w_c);
    }
}

// Template definition (a definition visible only in this TU)

template <bool ENCRYPT>

/**
 * @brief Efficiently applies the HPP rules for propagation, collision, and reflection using flat memory buffers.
 *
 * This function processes a single cell using fast pointer and stride arithmetic, suitable for high-performance parallel execution.
 * It supports both encryption (forward) and decryption (backward) modes, depending on the template parameter.
 * @param G Pointer to the flat memory buffer representing the grid.
 * @param rowStride The stride (width) of the local grid, including halos. (W = local_cols + 2)
 * @param i Local row index (1..local_rows), including halos.
 * @param j Local column index (1..local_cols), including halos.
 * @param wrow Pointer to the wall mask row for the current global row (column 0).
 * @param wrow_up Pointer to the wall mask row for the row above the current global row (row-1|mod Nmask).
 * @param wrow_dn Pointer to the wall mask row for the row below the current global row (row+1|mod Nmask).
 * @param Nmask The global width of the wall mask (grid_size).
 * @param j0_global The global column index that corresponds to local column j==1 (offset_cols).
 * @return The updated lattice after applying the HPP rules.
 */
inline uint8_t applyRules_fast(
    const uint8_t *FAST_RESTRICT G,
    int rowStride,
    int i, int j,
    const uint8_t *FAST_RESTRICT wrow,
    const uint8_t *FAST_RESTRICT wrow_up,
    const uint8_t *FAST_RESTRICT wrow_dn,
    int Nmask,
    int j0_global)
{
    // Local j-indexing (with halos) for the data buffer
    const int S = rowStride;                 // S = local_cols + 2
    const int jl = (j == 0 ? S - 1 : j - 1); // (becomes j=1 -> 0 = left halo)
    const int jr = (j == S - 1 ? 0 : j + 1); // (becomes j=local_cols -> local_cols+1 = right halo)

    // Line offsets in flat local memory
    const std::size_t row = static_cast<std::size_t>(i) * static_cast<std::size_t>(S);
    const std::size_t row_up = static_cast<std::size_t>(i - 1) * static_cast<std::size_t>(S);
    const std::size_t row_dn = static_cast<std::size_t>(i + 1) * static_cast<std::size_t>(S);

    // Read cells from the local buffer
    const uint8_t c = G[row + j];
    const uint8_t up = G[row_up + j];
    const uint8_t dn = G[row_dn + j];
    const uint8_t lf = G[row + jl];
    const uint8_t rt = G[row + jr];

    // --- Global column indexes for the wall mask ---
    // Local j (1..local_cols) corresponds to global: jg = j0_global + (j-1)
    // Halos j=0 / j=local_cols+1 -> jg-1 / jg+1 (with Wrap over Nmask)
    const int jg = (j0_global + (j - 1) + Nmask) % Nmask;
    const int jlg = (jg - 1 + Nmask) % Nmask;
    const int jrg = (jg + 1) % Nmask;

    // Wall flags (wrow* indicate the START of the respective global line, therefore direct indexing with jg/jlg/jrg)
    const bool w_c = (wrow[jg] != 0);
    const bool w_up = (wrow_up[jg] != 0);
    const bool w_dn = (wrow_dn[jg] != 0);
    const bool w_lf = (wrow[jlg] != 0);
    const bool w_rt = (wrow[jrg] != 0);

    if constexpr (ENCRYPT)
    {
        // ---------- FORWARD: collision -> propagate -> reflection ----------
        const uint8_t c_col = collision(c, w_c);
        const uint8_t up_col = collision(up, w_up);
        const uint8_t dn_col = collision(dn, w_dn);
        const uint8_t lf_col = collision(lf, w_lf);
        const uint8_t rt_col = collision(rt, w_rt);

        uint8_t next = static_cast<uint8_t>(c_col & 0xF0); // Keep high bits (upper 4 bits)
        if (up_col & 0x02)
            next |= 0x02; // SOUTH from ABOVE
        if (dn_col & 0x08)
            next |= 0x08; // NORTH from BELOW
        if (lf_col & 0x04)
            next |= 0x04; // EAST  from LEFT
        if (rt_col & 0x01)
            next |= 0x01; // WEST  from RIGHT

        return reflection(next, w_c);
    }
    else
    {
        // ---------- BACKWARD/INVERSE: inv_reflection -> inv_propagate -> inv_collision ----------
        const uint8_t c_ref = inverse_reflection(c, w_c);
        const uint8_t up_ref = inverse_reflection(up, w_up);
        const uint8_t dn_ref = inverse_reflection(dn, w_dn);
        const uint8_t lf_ref = inverse_reflection(lf, w_lf);
        const uint8_t rt_ref = inverse_reflection(rt, w_rt);

        uint8_t prev = static_cast<uint8_t>(c_ref & 0xF0);
        if (up_ref & 0x08)
            prev |= 0x08; // NORTH from ABOVE
        if (dn_ref & 0x02)
            prev |= 0x02; // SOUTH from BELOW
        if (rt_ref & 0x04)
            prev |= 0x04; // EAST  from RIGHT
        if (lf_ref & 0x01)
            prev |= 0x01; // WEST  from LEFT

        return inverse_collision(prev, w_c);
    }
}

// Explicit instantiations: exactly two versions (Encrypt & Decrypt)
template uint8_t applyRules_fast<true>(
    const uint8_t *FAST_RESTRICT, int, int, int,
    const uint8_t *FAST_RESTRICT, const uint8_t *FAST_RESTRICT, const uint8_t *FAST_RESTRICT,
    int, int);
template uint8_t applyRules_fast<false>(
    const uint8_t *FAST_RESTRICT, int, int, int,
    const uint8_t *FAST_RESTRICT, const uint8_t *FAST_RESTRICT, const uint8_t *FAST_RESTRICT,
    int, int);
