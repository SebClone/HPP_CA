// Implementiert die Funktionen, die in hpp_encryptor.hpp festgelegt sind
#include "hpp_encryptor.hpp"

#include <fstream>
#include <random>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>

// Für die Bildverarbeitung
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

void printBits(uint8_t value)
{
    std::bitset<8> bits(value);
    std::cout << bits;
}

void printGrid(const std::vector<std::vector<uint8_t>> &grid)
{
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

std::vector<uint8_t> readFileBytes(const std::string &filename)
{
    std::ifstream file(filename, std::ios::binary); // open file in binary mode
    std::vector<uint8_t> data;

    if (!file)
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    // Move to the end to get file size
    file.seekg(0, std::ios::end);
    size_t size = file.tellg();

    // Resize vector to fit file contents
    data.resize(size);

    // Go back to start and read data into vector
    file.seekg(0);
    file.read(reinterpret_cast<char *>(data.data()), size);

    return data;
}

// Lädt ein Bild belibiges .png/ .jpg Bild als RGB-Daten
std::vector<uint8_t> loadImageAsRGB(const std::string &filename, int &width, int &height, int &channels)
{
    unsigned char *data = stbi_load(filename.c_str(), &width, &height, &channels, 3);
    if (!data)
    {
        std::cerr << "Fehler: Bild konnte nicht geladen werden: " << filename << std::endl;
        return {};
    }

    size_t size = width * height * 3; // 3 Bytes pro Pixel (R, G, B)
    std::vector<uint8_t> imgData(data, data + size);

    stbi_image_free(data); // Speicher freigeben
    return imgData;
}

std::vector<std::vector<uint8_t>> reshapeToMatrix(const std::vector<uint8_t> &data, size_t &grid_size)
{
    size_t originalSize = data.size();

    // Calculate the smallest square dimension that can hold all the data
    grid_size = static_cast<size_t>(std::ceil(std::sqrt(originalSize)));
    size_t paddedSize = grid_size * grid_size;

    // Copy and pad the data with 0s if needed
    std::vector<uint8_t> padded = data;
    padded.resize(paddedSize, 0);
    // Fill a 2D matrix with the padded data
    std::vector<std::vector<uint8_t>> matrix(grid_size, std::vector<uint8_t>(grid_size));
    for (size_t i = 0; i < paddedSize; ++i)
    {
        matrix[i / grid_size][i % grid_size] = padded[i];
    }

    return matrix;
}

std::vector<uint8_t> flattenMatrix(const std::vector<std::vector<uint8_t>> &matrix)
{
    std::vector<uint8_t> flat;
    for (const auto &row : matrix)
    {
        for (uint8_t value : row)
        {
            flat.push_back(value);
        }
    }
    return flat;
}

void saveAsAsciiText(const std::vector<uint8_t> &data, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (uint8_t byte : data)
    {
        file << static_cast<char>(byte); // Convert uint8_t to ASCII character
    }

    file.close();
}

// Speichert einen kompletten Zellautomaten-Frame als Binärdatei.
void save_frame_bin(const std::vector<uint8_t> &frame, int iter)
{
    // String-Stream, um den Dateinamen im gewünschten Format zu bauen
    std::ostringstream oss;
    oss << "frames/frame_" << std::setw(6) << std::setfill('0') << iter << ".bin";

    std::string filename = oss.str();
    // Öffne die Datei im Binärmodus
    std::ofstream ofs(filename, std::ios::binary);

    // Falls das Öffnen fehlschlägt, Fehlermeldung ausgeben und abbrechen
    if (!ofs)
    {
        std::cerr << "ERROR: could not open file '" << filename << "' for writing.\n";
        return;
    }

    // Schreibe den kompletten Frame als Bytes in die Datei
    ofs.write(reinterpret_cast<const char *>(frame.data()), frame.size());

    // Prüfe, ob beim Schreiben ein Fehler auftrat
    if (!ofs)
    {
        std::cerr << "ERROR: failed to write frame to '" << filename << "'\n";
    }
}

Mask generateRandomWallMask(int grid_size, double wall_ratio, uint32_t seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    Mask wall_mask(grid_size, std::vector<uint8_t>(grid_size, 0));
    int wall_count = 0;

    for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
        {
            if (dist(rng) < wall_ratio)
            {
                wall_mask[i][j] = 1;
                wall_count++;
            }
        }
    }

    std::cout << "Wall cells placed: " << wall_count << std::endl;
    return wall_mask;
}

void saveWallMaskBinary(const Mask &wall_mask, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Cannot write wall mask to " << filename << std::endl;
        return;
    }

    for (const auto &row : wall_mask)
    {
        for (uint8_t wall : row)
        {
            file.write(reinterpret_cast<char *>(&wall), 1);
        }
    }
}

Mask loadWallMaskBinary(int grid_size, const std::string &filename)
{
    Mask wall_mask(grid_size, std::vector<uint8_t>(grid_size, 0));
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Cannot read wall mask from " << filename << std::endl;
        return wall_mask;
    }

    for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
        {
            char wall;
            file.read(&wall, 1);
            wall_mask[i][j] = (wall == 1);
        }
    }

    return wall_mask;
}

// ----------------------------------------------
// Forward declarations of functions
// ----------------------------------------------
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

// ----------------------------------------------
// Backward declarations of functions
// ----------------------------------------------
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

/*
 *-------------------------------------------------------
 * Functions for MPI communication
 * -------------------------------------------------------
 */

void broadcastMask(Mask &mask, MPI_Comm comm)
{
    int rows = mask.size();
    int cols = rows > 0 ? mask[0].size() : 0;
    MPI_Bcast(&rows, 1, MPI_INT, 0, comm);
    MPI_Bcast(&cols, 1, MPI_INT, 0, comm);

    // Auf Nicht-Root: Speicher für mask anlegen
    mask.resize(rows);
    for (int i = 0; i < rows; ++i)
    {
        mask[i].resize(cols);
    }

    for (int i = 0; i < rows; ++i)
    {
        MPI_Bcast(mask[i].data(), cols, MPI_BYTE, 0, comm);
    }
}

void copy_n_bytes(const uint8_t *src, std::size_t count, uint8_t *dst)
{
    for (std::size_t i = 0; i < count; ++i)
    {
        dst[i] = src[i];
    }
}

// Gekapselte Hauptlogik der HPP-Operationen
uint8_t applyRules(
    const Matrix &active_grid,
    const Mask &wall_mask,
    bool doEncrypt,  // true = vorwärts (Encrypt), false = rückwärts (Decrypt)
    int i,           // lokale Zeile im Subgitter (mit Halo): 1..local_rows
    int j,           // Spalte: 0..N-1
    int offset_rows) // globale Startzeile dieses Ranks (ohne Halo)
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
    const uint8_t c = active_grid[i][j];
    const uint8_t up = active_grid[iu][j];
    const uint8_t dn = active_grid[id][j];
    const uint8_t lf = active_grid[i][jl];
    const uint8_t rt = active_grid[i][jr];

    // Wand-Flags (Maske ist uint8_t → nonzero == true)
    const bool w_c = wall_mask[gr][gc] != 0;
    const bool w_up = wall_mask[(gr - 1 + N) % N][gc] != 0;
    const bool w_dn = wall_mask[(gr + 1) % N][gc] != 0;
    const bool w_lf = wall_mask[gr][jl] != 0;
    const bool w_rt = wall_mask[gr][jr] != 0;

    if (doEncrypt)
    {
        // ---------- VORWÄRTS: collision -> propagate (in diese Zielzelle sammeln) -> reflection ----------

        // 1) collision auf Zentrum + Nachbarn
        const uint8_t c_col = collision(c, w_c);
        const uint8_t up_col = collision(up, w_up);
        const uint8_t dn_col = collision(dn, w_dn);
        const uint8_t lf_col = collision(lf, w_lf);
        const uint8_t rt_col = collision(rt, w_rt);

        // 2) "propagation" in diese Zelle:
        //    obere 4 Bits vom Zentrum behalten, Richtungsbits aus kollidierten Nachbarn einsammeln
        uint8_t next = static_cast<uint8_t>(c_col & 0b11110000);

        if (up_col & 0b00000010)
            next |= 0b00000010; // von oben kommt dessen SOUTH
        if (dn_col & 0b00001000)
            next |= 0b00001000; // von unten kommt dessen NORTH
        if (lf_col & 0b00000100)
            next |= 0b00000100; // von links kommt dessen EAST
        if (rt_col & 0b00000001)
            next |= 0b00000001; // von rechts kommt dessen WEST

        // 3) reflection auf der Zielzelle
        return reflection(next, w_c);
    }
    else
    {
        // ---------- RÜCKWÄRTS: inverse_reflection -> inverse_propagate (sammeln) -> inverse_collision ----------

        // 1) inverse_reflection auf Zentrum + Nachbarn
        const uint8_t c_ref = inverse_reflection(c, w_c);
        const uint8_t up_ref = inverse_reflection(up, w_up);
        const uint8_t dn_ref = inverse_reflection(dn, w_dn);
        const uint8_t lf_ref = inverse_reflection(lf, w_lf);
        const uint8_t rt_ref = inverse_reflection(rt, w_rt);

        // 2) inverse_propagate „in diese Zelle“ (entspricht deiner grid-weiten Variante):
        uint8_t prev = static_cast<uint8_t>(c_ref & 0b11110000);

        if (up_ref & 0b00001000)
            prev |= 0b00001000; // N  kommt von OBEN
        if (dn_ref & 0b00000010)
            prev |= 0b00000010; // S  kommt von UNTEN
        if (rt_ref & 0b00000100)
            prev |= 0b00000100; // E  kommt von RECHTS
        if (lf_ref & 0b00000001)
            prev |= 0b00000001; // W  kommt von LINKS

        // 3) inverse_collision
        return inverse_collision(prev, w_c);
    }
}

void saveBinary(const std::vector<uint8_t> &data, const char *filename)
{
    std::ofstream out(filename, std::ios::binary);
    out.write(reinterpret_cast<const char *>(data.data()), static_cast<std::streamsize>(data.size()));
}

void saveEncryptedMeta(uint64_t originalSize, uint32_t gridSize, const char *metaFilename)
{
    std::ofstream out(metaFilename, std::ios::binary);
    const uint32_t magic = 0x48505031; // "HPP1"
    const uint32_t version = 1;
    out.write(reinterpret_cast<const char *>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char *>(&version), sizeof(version));
    out.write(reinterpret_cast<const char *>(&originalSize), sizeof(originalSize));
    out.write(reinterpret_cast<const char *>(&gridSize), sizeof(gridSize));
}

bool loadEncryptedMeta(uint64_t &originalSize, uint32_t &gridSize, const char *metaFilename)
{
    std::ifstream in(metaFilename, std::ios::binary);
    if (!in)
        return false;
    uint32_t magic = 0, version = 0;
    in.read(reinterpret_cast<char *>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char *>(&version), sizeof(version));
    if (magic != 0x48505031 || version != 1)
        return false;
    in.read(reinterpret_cast<char *>(&originalSize), sizeof(originalSize));
    in.read(reinterpret_cast<char *>(&gridSize), sizeof(gridSize));
    return static_cast<bool>(in);
}