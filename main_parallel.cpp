#include "hpp_encryptor.hpp"
#include <cmath>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdio>
#include <deque>      // für Pending-Requests
#include <algorithm>  // std::fill, std::min

using Matrix = std::vector<std::vector<uint8_t>>;

constexpr int NUM_ITERATIONS = 1000;
constexpr int TAG_NS = 0;

static constexpr const char* ENC_BIN  = "encrypted_full.bin";
static constexpr const char* ENC_META = "encrypted_full.meta";
static constexpr const char* PLAIN_IN = "message.txt";
static constexpr const char* DECRYPT_OUT = "decrypted_message.txt";

#ifndef DUMP_FRAMES
    #define DUMP_FRAMES 1
#endif

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Modus: true = Verschlüsseln, false = Entschlüsseln
    bool doEncrypt = false;

    int mode = doEncrypt ? 1 : 0;
    MPI_Bcast(&mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    doEncrypt = (mode != 0);

    const char *key_filename = "wall_mask.key";

    double t_start = MPI_Wtime();

    std::vector<uint8_t> fileData;
    uint64_t originalSize = 0;
    int grid_size = 0;

    // --- PORTABILITY FIX: finde einen 64-bit Integer-MPI-Datentyp ---
    MPI_Datatype MPI_UINT64_MATCHED;
    {
        int rc = MPI_Type_match_size(MPI_TYPECLASS_INTEGER, 8, &MPI_UINT64_MATCHED);
        if (rc != MPI_SUCCESS) {
            // Fallback – sollte praktisch nie passieren; prüfe Größe trotzdem
            static_assert(sizeof(unsigned long long) == 8, "Need 64-bit unsigned long long");
            MPI_UINT64_MATCHED = MPI_UNSIGNED_LONG_LONG;
        }
    }

    if (rank == 0) {
        if (doEncrypt) {
            MPI_File fhin;
            if (MPI_File_open(MPI_COMM_SELF, PLAIN_IN, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhin) != MPI_SUCCESS) {
                std::cerr << "ERROR: could not open input file '" << PLAIN_IN << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            MPI_Offset fsz = 0;
            MPI_File_get_size(fhin, &fsz);
            MPI_File_close(&fhin);

            if (fsz < 0) {
                std::cerr << "ERROR: invalid size for '" << PLAIN_IN << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            originalSize = static_cast<uint64_t>(fsz);
            grid_size = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(originalSize))));
            std::cout << "Encrypting mode selected.\n";
            std::cout << "Read size " << originalSize << " bytes from " << PLAIN_IN
                      << ", grid size = " << grid_size << std::endl;
        } 
        else {
            // Meta laden (liefert Originalgröße & N)
            uint32_t Nmeta = 0;
            if (!loadEncryptedMeta(originalSize, Nmeta, ENC_META)) {
                std::cerr << "ERROR: could not read meta file '" << ENC_META << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            grid_size = static_cast<int>(Nmeta);

            std::cout << "Decrypting mode selected.\n";
            std::cout << "Loaded meta from " << ENC_META << ", grid size = " << grid_size
                      << ", originalSize = " << originalSize << std::endl;
        }
    }

    MPI_Bcast(&originalSize, 1, MPI_UINT64_MATCHED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&grid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

     if (nprocs > grid_size) {
        if (rank == 0) std::cerr << "nprocs > grid_size: mindestens ein Rank bekäme 0 Zeilen – für Torus ungültig.\n";
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    int rows_per_rank = grid_size / nprocs;
    int remainder = grid_size % nprocs;
    int local_rows = rows_per_rank + (rank < remainder ? 1 : 0);
    int offset_rows = rank * rows_per_rank + std::min(rank, remainder);

    size_t paddedSize = static_cast<size_t>(grid_size) * grid_size;
    std::vector<uint8_t> local_core(local_rows * grid_size);

    if (doEncrypt) {
        MPI_File fh;
        if (MPI_File_open(MPI_COMM_WORLD, PLAIN_IN, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
            if (rank == 0) std::cerr << "ERROR: could not open input file '" << PLAIN_IN << "'\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        MPI_File_set_atomicity(fh, 0);

        // Puffer initial auf 0 (Padding lokal)
        std::fill(local_core.begin(), local_core.end(), 0);

        const MPI_Offset base = static_cast<MPI_Offset>(offset_rows) * grid_size;
        const size_t want = static_cast<size_t>(local_rows) * grid_size;

        size_t remaining = 0;
        if (static_cast<uint64_t>(base) < originalSize) {
            remaining = static_cast<size_t>(
                std::min<uint64_t>(want, originalSize - static_cast<uint64_t>(base))
            );
        }

        if (remaining > 0) {
            MPI_File_read_at_all(fh, base, local_core.data(),
                                 static_cast<int>(remaining), MPI_BYTE, MPI_STATUS_IGNORE);
        }
        MPI_File_close(&fh);
    }
    else {
        // Jeder Rank liest seinen zusammenhängenden Zeilenblock aus ENC_BIN.
        MPI_File fh;
        if (MPI_File_open(MPI_COMM_WORLD, ENC_BIN, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
            if (rank == 0) std::cerr << "ERROR: could not open '" << ENC_BIN << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Optional: Atomicity aus (typischerweise schneller)
        MPI_File_set_atomicity(fh, 0);

         MPI_Offset fsz = 0;
        MPI_File_get_size(fh, &fsz);
        if (fsz != static_cast<MPI_Offset>(paddedSize)) {
            if (rank == 0) {
                std::cerr << "ERROR: '" << ENC_BIN << "' has " << fsz
                          << " bytes, expected " << paddedSize << "\n";
            }
            MPI_File_close(&fh);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        MPI_Offset file_off = static_cast<MPI_Offset>(offset_rows) * grid_size; // Byte-Offset, da 1 Byte/Element
        MPI_File_read_at_all(fh, file_off,
                             local_core.data(),
                             static_cast<int>(local_core.size()),
                             MPI_BYTE, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
    }
    
    Matrix grid_current(local_rows + 2, std::vector<uint8_t>(grid_size));
    Matrix grid_next = grid_current;

    for (int i = 0; i < local_rows; ++i) {
        copy_n_bytes(
            local_core.data() + i * grid_size,
            static_cast<std::size_t>(grid_size),
            grid_current[i + 1].data());
    }

    Mask wall_mask;
    if (rank == 0) {
        if (doEncrypt) {
            wall_mask = generateRandomWallMask(grid_size, 0.1);
            saveWallMaskBinary(wall_mask, key_filename);
            std::cout << "Wall mask generated and saved to " << key_filename << "\n";
        } 
        else {
            wall_mask = loadWallMaskBinary(grid_size, key_filename);
            std::cout << "Wall mask loaded from " << key_filename << "\n";
        }
    }
    broadcastMask(wall_mask, MPI_COMM_WORLD);


    // MPI-Topologie für Nachbarschafts-Kommunikation
    int dims[1]     = { nprocs }; // 1D-Topologie mit nprocs Prozessen
    int periods[1]  = { 1 };      // 1 = periodisch (Torus), 0 = nicht periodisch
    int reorder     = 0;          // MPI darf Ranks nicht neu anordnen für bessere Nachbarschaft  
    MPI_Comm cart_comm; 
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &cart_comm);

    int up, down;
    MPI_Cart_shift(cart_comm, 0, 1, &up, &down);

    MPI_Request halo_reqs_current[4], halo_reqs_next[4];

    MPI_Recv_init(grid_current[0].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_current[0]);
    MPI_Recv_init(grid_current[local_rows + 1].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_current[1]);
    MPI_Send_init(grid_current[1].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_current[2]);
    MPI_Send_init(grid_current[local_rows].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_current[3]);

    MPI_Recv_init(grid_next[0].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_next[0]);
    MPI_Recv_init(grid_next[local_rows + 1].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_next[1]);
    MPI_Send_init(grid_next[1].data(), grid_size, MPI_UINT8_T, up, TAG_NS, cart_comm, &halo_reqs_next[2]);
    MPI_Send_init(grid_next[local_rows].data(), grid_size, MPI_UINT8_T, down, TAG_NS, cart_comm, &halo_reqs_next[3]);

    Matrix* active_grid = &grid_current;
    Matrix* target_grid = &grid_next;
    MPI_Request* active_reqs = halo_reqs_current;
    MPI_Request* inactive_reqs = halo_reqs_next;

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        MPI_Startall(4, active_reqs);

        // Innenbereich berechnen
        if (local_rows >= 3) {
            for (int i = 2; i <= local_rows - 1; i++) {
                for (int j = 0; j < grid_size; ++j) {
                    (*target_grid)[i][j] = applyRules(*active_grid, wall_mask, doEncrypt, i, j, offset_rows);
                }
            }
        }

        MPI_Waitall(4, active_reqs, MPI_STATUSES_IGNORE);

        if (local_rows >= 1) {
            for (int j = 0; j < grid_size; ++j) {
                (*target_grid)[1][j] = applyRules(*active_grid, wall_mask, doEncrypt, 1, j, offset_rows);
            }
        }
        if (local_rows >= 2) {
            for (int j = 0; j < grid_size; ++j) {
                (*target_grid)[local_rows][j] = applyRules(*active_grid, wall_mask, doEncrypt, local_rows, j, offset_rows);
            }
        }
        std::swap(active_grid, target_grid);
        std::swap(active_reqs, inactive_reqs);

        // Frames speichern
        #if DUMP_FRAMES
        if (iter % 5 == 0 || iter == NUM_ITERATIONS - 1) {
            // lokalen Frame-Puffer bilden
            std::vector<uint8_t> frame_core(local_rows * grid_size);
            for (int i = 0; i < local_rows; ++i) {
                copy_n_bytes(
                    (*active_grid)[i + 1].data(),
                    static_cast<std::size_t>(grid_size),
                    frame_core.data() + i * grid_size);
            }

            // Eine Datei pro Frame
            char fname[128];
            std::snprintf(fname, sizeof(fname), "frame_%05d.bin", iter);

            MPI_File fhf;
            MPI_File_open(MPI_COMM_WORLD, fname,
                          MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fhf);
            MPI_File_set_atomicity(fhf, 0);

            // Jeder Rank schreibt an seinen Zeilen-Offset
            MPI_Offset off = static_cast<MPI_Offset>(offset_rows) * grid_size;

            // Nicht-blockierendes kollektives I/O (kann man ohne sofortiges Wait überlappen)
            MPI_Request req;
            MPI_File_iwrite_at_all(fhf, off,
                                   frame_core.data(),
                                   static_cast<int>(frame_core.size()),
                                   MPI_BYTE, &req);

            // Für Einfachheit warten wir direkt; für Overlap: Requests sammeln und später warten
            MPI_Wait(&req, MPI_STATUS_IGNORE);
            MPI_File_close(&fhf);
        }
        #endif
    }
    for (auto &r : halo_reqs_current) MPI_Request_free(&r);
    for (auto &r : halo_reqs_next)    MPI_Request_free(&r);
    MPI_Comm_free(&cart_comm);

    std::vector<uint8_t> result_core(local_rows * grid_size);
    for (int i = 0; i < local_rows; ++i)
    {
        copy_n_bytes(
            (*active_grid)[i + 1].data(),
            static_cast<std::size_t>(grid_size),
            result_core.data() + i * grid_size);
    }

    if (doEncrypt) {  
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, ENC_BIN,
                    MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_atomicity(fh, 0);
        if (rank == 0) {
            MPI_File_set_size(fh, static_cast<MPI_Offset>(paddedSize));
        }
        MPI_Barrier(MPI_COMM_WORLD); 

        MPI_Offset file_off = static_cast<MPI_Offset>(offset_rows) * grid_size;
        MPI_File_write_at_all(fh, file_off,
                              result_core.data(),
                              static_cast<int>(result_core.size()),
                              MPI_BYTE, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);

        if (rank == 0) {
            // Meta bleibt klein → auf Rank 0 schreiben
            saveEncryptedMeta(originalSize, static_cast<uint32_t>(grid_size), ENC_META);
            std::cout << "Saved " << ENC_BIN << " (parallel, "
                      << paddedSize << " bytes) and " << ENC_META
                      << " (originalSize=" << originalSize
                      << ", grid_size=" << grid_size << ")\n";
        }
    } 
    else {
        // *** CHANGE: Decrypt -> Klartext parallel schreiben (ohne Gather)
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, DECRYPT_OUT,
                      MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_atomicity(fh, 0);

        if (rank == 0) {
            // Datei genau auf originalSize einstellen (kein Padding im Output)
            MPI_File_set_size(fh, static_cast<MPI_Offset>(originalSize));
        }
        MPI_Barrier(MPI_COMM_WORLD);

        const MPI_Offset base = static_cast<MPI_Offset>(offset_rows) * grid_size;
        const size_t have = static_cast<size_t>(local_rows) * grid_size;

        size_t to_write = 0;
        if (static_cast<uint64_t>(base) < originalSize) {
            to_write = static_cast<size_t>(
                std::min<uint64_t>(have, originalSize - static_cast<uint64_t>(base))
            );
        }

        if (to_write > 0) {
            MPI_File_write_at_all(fh, base, result_core.data(),
                                  static_cast<int>(to_write), MPI_BYTE, MPI_STATUS_IGNORE);
        }
        MPI_File_close(&fh);

        if (rank == 0) {
            std::cout << "Saved " << DECRYPT_OUT << " (" << originalSize << " bytes, parallel)\n";
        }   
    }

    double t_end = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << std::fixed << std::setprecision(6)
                  << "Total runtime: " << (t_end - t_start) << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}

