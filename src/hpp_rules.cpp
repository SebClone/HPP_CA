#include "hpp_rules.hpp"

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

uint8_t inverse_reflection(uint8_t current_cell, bool is_wall)
{
    return reflection(current_cell, is_wall);
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

uint8_t inverse_collision(uint8_t current_cell, bool is_wall)
{
    return collision(current_cell, is_wall);
}


// Gekapselte Hauptlogik der HPP-Operationen
uint8_t applyRules(
    const Matrix& active_grid,
    const Mask&   wall_mask,
    bool          doEncrypt,   // true = vorwärts (Encrypt), false = rückwärts (Decrypt)
    int           i,           // lokale Zeile im Subgitter (mit Halo): 1..local_rows
    int           j,           // Spalte: 0..N-1
    int           offset_rows) // globale Startzeile dieses Ranks (ohne Halo)
{
    const int N = static_cast<int>(wall_mask.size()); // quadratisches Gitter

    // Globale Koordinate der Zielzelle (Zeilen über Rank, Spalten lokal im Torus)
    const int gr = (offset_rows + (i - 1) + N) % N;
    const int gc = j;

    // Nachbarindizes (lokal & global)
    const int iu = i - 1;           // Nachbar oben   (kommt aus Halo/Innen)
    const int id = i + 1;           // Nachbar unten
    const int jl = (j - 1 + N) % N; // links  (torus)
    const int jr = (j + 1) % N;     // rechts (torus)

    // Zellen laden (nur lesen!)
    const uint8_t c  = active_grid[i ][j ];
    const uint8_t up = active_grid[iu][j ];
    const uint8_t dn = active_grid[id][j ];
    const uint8_t lf = active_grid[i ][jl];
    const uint8_t rt = active_grid[i ][jr];

    // Wand-Flags (Maske ist uint8_t → nonzero == true)
    const bool w_c  = wall_mask[gr][gc] != 0;
    const bool w_up = wall_mask[(gr - 1 + N) % N][gc] != 0;
    const bool w_dn = wall_mask[(gr + 1) % N][gc]     != 0;
    const bool w_lf = wall_mask[gr][jl]               != 0;
    const bool w_rt = wall_mask[gr][jr]               != 0;

    if (doEncrypt)
    {
        // ---------- VORWÄRTS: collision -> propagate (in diese Zielzelle sammeln) -> reflection ----------

        // 1) collision auf Zentrum + Nachbarn
        const uint8_t c_col  = collision(c , w_c );
        const uint8_t up_col = collision(up, w_up);
        const uint8_t dn_col = collision(dn, w_dn);
        const uint8_t lf_col = collision(lf, w_lf);
        const uint8_t rt_col = collision(rt, w_rt);

        // 2) "propagation" in diese Zelle:
        //    obere 4 Bits vom Zentrum behalten, Richtungsbits aus kollidierten Nachbarn einsammeln
        uint8_t next = static_cast<uint8_t>(c_col & 0b11110000);

        if (up_col & 0b00000010) next |= 0b00000010; // von oben kommt dessen SOUTH
        if (dn_col & 0b00001000) next |= 0b00001000; // von unten kommt dessen NORTH
        if (lf_col & 0b00000100) next |= 0b00000100; // von links kommt dessen EAST
        if (rt_col & 0b00000001) next |= 0b00000001; // von rechts kommt dessen WEST

        // 3) reflection auf der Zielzelle
        return reflection(next, w_c);
    }
    else
    {
        // ---------- RÜCKWÄRTS: inverse_reflection -> inverse_propagate (sammeln) -> inverse_collision ----------

        // 1) inverse_reflection auf Zentrum + Nachbarn
        const uint8_t c_ref  = inverse_reflection(c , w_c );
        const uint8_t up_ref = inverse_reflection(up, w_up);
        const uint8_t dn_ref = inverse_reflection(dn, w_dn);
        const uint8_t lf_ref = inverse_reflection(lf, w_lf);
        const uint8_t rt_ref = inverse_reflection(rt, w_rt);

        // 2) inverse_propagate in Zelle (entspricht grid-weiten Variante):
        uint8_t prev = static_cast<uint8_t>(c_ref & 0b11110000);

        if (up_ref & 0b00001000) prev |= 0b00001000; // N  kommt von OBEN
        if (dn_ref & 0b00000010) prev |= 0b00000010; // S  kommt von UNTEN
        if (rt_ref & 0b00000100) prev |= 0b00000100; // E  kommt von RECHTS
        if (lf_ref & 0b00000001) prev |= 0b00000001; // W  kommt von LINKS

        // 3) inverse_collision
        return inverse_collision(prev, w_c);
    }
}



// Template-Definition (eine Definition, sichtbar nur in dieser TU)
template<bool ENCRYPT>
inline uint8_t applyRules_fast(
    const uint8_t* FAST_RESTRICT G,
    int N,
    int i, int j,
    const uint8_t* FAST_RESTRICT wrow,
    const uint8_t* FAST_RESTRICT wrow_up,
    const uint8_t* FAST_RESTRICT wrow_dn
)
{
    // Toroidales j-Indexing (ohne % für j±1)
    const int jl = (j == 0   ? N - 1 : j - 1);
    const int jr = (j == N-1 ? 0     : j + 1);

    // Zeilenoffsets im flachen Speicher
    const std::size_t row    = static_cast<std::size_t>(i)     * static_cast<std::size_t>(N);
    const std::size_t row_up = static_cast<std::size_t>(i - 1) * static_cast<std::size_t>(N);
    const std::size_t row_dn = static_cast<std::size_t>(i + 1) * static_cast<std::size_t>(N);

    // Zellen lesen
    const uint8_t c  = G[row    + j];
    const uint8_t up = G[row_up + j];
    const uint8_t dn = G[row_dn + j];
    const uint8_t lf = G[row    + jl];
    const uint8_t rt = G[row    + jr];

    // Wand-Flags (Zeiger auf fertige Wall-Zeilen kommen von außen)
    const bool w_c  = (wrow   [j ] != 0);
    const bool w_up = (wrow_up[j ] != 0);
    const bool w_dn = (wrow_dn[j ] != 0);
    const bool w_lf = (wrow   [jl] != 0);
    const bool w_rt = (wrow   [jr] != 0);

    if constexpr (ENCRYPT) {
        // ---------- vorwärts: collision -> propagate -> reflection ----------
        const uint8_t c_col  = collision(c , w_c );
        const uint8_t up_col = collision(up, w_up);
        const uint8_t dn_col = collision(dn, w_dn);
        const uint8_t lf_col = collision(lf, w_lf);
        const uint8_t rt_col = collision(rt, w_rt);

        uint8_t next = static_cast<uint8_t>(c_col & 0xF0); // obere 4 Bits behalten
        if (up_col & 0x02) next |= 0x02; // SOUTH von oben
        if (dn_col & 0x08) next |= 0x08; // NORTH von unten
        if (lf_col & 0x04) next |= 0x04; // EAST  von links
        if (rt_col & 0x01) next |= 0x01; // WEST  von rechts

        return reflection(next, w_c);
    } else {
        // ---------- rückwärts: inv_reflection -> inv_propagate -> inv_collision ----------
        const uint8_t c_ref  = inverse_reflection(c , w_c );
        const uint8_t up_ref = inverse_reflection(up, w_up);
        const uint8_t dn_ref = inverse_reflection(dn, w_dn);
        const uint8_t lf_ref = inverse_reflection(lf, w_lf);
        const uint8_t rt_ref = inverse_reflection(rt, w_rt);

        uint8_t prev = static_cast<uint8_t>(c_ref & 0xF0);
        if (up_ref & 0x08) prev |= 0x08; // N  von oben
        if (dn_ref & 0x02) prev |= 0x02; // S  von unten
        if (rt_ref & 0x04) prev |= 0x04; // E  von rechts
        if (lf_ref & 0x01) prev |= 0x01; // W  von links

        return inverse_collision(prev, w_c);
    }
}

// Explizite Instanziierungen: erzeugen genau zwei Versionen (Encrypt & Decrypt)
// Dadurch brauchen andere Übersetzungseinheiten nur den Header – sie instanziieren NICHT erneut.
template uint8_t applyRules_fast<true>(
    const uint8_t* FAST_RESTRICT, int, int, int,
    const uint8_t* FAST_RESTRICT, const uint8_t* FAST_RESTRICT, const uint8_t* FAST_RESTRICT
);
template uint8_t applyRules_fast<false>(
    const uint8_t* FAST_RESTRICT, int, int, int,
    const uint8_t* FAST_RESTRICT, const uint8_t* FAST_RESTRICT, const uint8_t* FAST_RESTRICT
);