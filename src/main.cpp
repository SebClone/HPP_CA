#include "hpp_rules.hpp"
#include "app_config.hpp"
#include "io.hpp"
#include "utilities.hpp"

#include <mpi.h>        // MPI_Init, MPI_Comm_*, MPI_* APIs
#include <vector>       // std::vector
#include <string>       // std::string
#include <iostream>     // std::cout, std::cerr
#include <iomanip>      // std::fixed, std::setprecision
#include <algorithm>    // std::min
#include <cmath>        // std::ceil, std::sqrt
#include <cstdio>       // std::snprintf
#include <cstdint>      // std::uint64_t, std::uint8_t
#include <cstddef>      // std::size_t
#include <cstring>      // std::memcpy

using Matrix = std::vector<std::vector<uint8_t>>;
using Mask = std::vector<std::vector<uint8_t>>;

constexpr int TAG_NS = 0;

int main(int argc, char **argv) {

    // Initialisieren von MPI und Config
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    AppConfig cfg = parse_args_or_default(argc, argv, rank);
    bcast_config(cfg, /*root=*/0, MPI_COMM_WORLD);

    // Configs: Verhalten
    const bool doEncrypt      = (cfg.mode == AppMode::Encrypt); // true = Encrypt, false = Decrypt
    const int numIterations   = cfg.iterations;                 // Anzahl Iterationen
    const int frameInterval   = cfg.frame_interval;             // Frames speichern alle X Iterationen
    const bool dumpFrames     = cfg.dump_frames;                // true = Frames speichern, false = Frames nicht speichern
    const bool cartReorder    = cfg.cart_reorder;               // true = MPI_Cart_create darf Ranks neu anordnen
    const bool atomicIO       = cfg.atomic_io;                  // true = MPI_File_set_atomicity, false = nicht atomar
    const double wallDensity  = cfg.wall_density;               // % Wandzellen (0.0 - 1.0)
    const std::uint64_t seed  = cfg.seed;                       // seed für Wandzellen-Generierung (0 = zufälliger seed)

    // Configs: Eingabe- und Ausgabe
    const std::string& inPlain   = cfg.input;
    const std::string& encBin    = cfg.enc_bin;
    const std::string& metaPath  = cfg.meta;
    const std::string& keyPath   = cfg.key;
    const std::string& outPlain  = cfg.output;

    // Ausgeben der Konfiguration
    if (rank == 0) {
        std::cout << "Configuration:\n";
        std::cout << "Mode: " << (doEncrypt ? "Encrypt" : "Decrypt") << "\n";
        std::cout << "Iterations: " << numIterations << "\n";
        std::cout << "Grid Size: " << cfg.grid_size << "\n";
        std::cout << "Frame Interval: " << frameInterval << "\n";
        std::cout << "Dump Frames: " << (dumpFrames ? "Yes" : "No") << "\n";
        std::cout << "Atomic I/O: " << (atomicIO ? "Yes" : "No") << "\n";
        std::cout << "Cart Reorder: " << (cartReorder ? "Yes" : "No") << "\n";
        std::cout << "Wall Density: " << wallDensity * 100.0 << "%\n";
        std::cout << "Seed: " << seed << "\n";
        std::cout << "Input File: " << inPlain << "\n";
        std::cout << "Encrypted Binary File: " << encBin << "\n";
        std::cout << "Meta File: " << metaPath << "\n";
        std::cout << "Key File: " << keyPath << "\n";
        std::cout << "Output File: " << outPlain << "\n";
    }

    double t_start = MPI_Wtime();

    uint64_t originalSize = 0;
    int grid_size = 0;



    // Findet einen zu 64 Bit passenden MPI-Integer-Typ (für Portierbarkeit der originalSize Variablen)
    MPI_Datatype MPI_UINT64_MATCHED;
    int rc = MPI_Type_match_size(MPI_TYPECLASS_INTEGER, 8, &MPI_UINT64_MATCHED);
    if (rc != MPI_SUCCESS) {
        static_assert(sizeof(unsigned long long) == 8, "Need 64-bit unsigned long long");
        MPI_UINT64_MATCHED = MPI_UNSIGNED_LONG_LONG;
    }

    // Bestimmen der grid_size und der originalSize
    if (rank == 0) {
        if (doEncrypt) {
            // Encrypt: Rank 0 liest die Eingabedatei und bestimmt die Größe
            MPI_File fhin;
            if (MPI_File_open(MPI_COMM_SELF, inPlain.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fhin) != MPI_SUCCESS) {
                std::cerr << "ERROR: could not open input file '" << inPlain.c_str() << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            MPI_Offset fsz = 0;
            MPI_File_get_size(fhin, &fsz);
            MPI_File_close(&fhin);

            if (fsz < 0) {
                std::cerr << "ERROR: invalid size for '" << inPlain.c_str() << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            originalSize = static_cast<uint64_t>(fsz);
            grid_size = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(originalSize))));
            std::cout << "Encrypting mode selected.\n";
            std::cout << "Read size " << originalSize << " bytes from " << inPlain.c_str()
                      << ", grid size = " << grid_size << std::endl;
        } 
        else {
            // Decrypt: Rank 0 liest die Meta-Datei und bestimmt die Größe
            uint32_t Nmeta = 0;
            if (!read_meta_rank0(metaPath, originalSize, Nmeta)) {
                std::cerr << "ERROR: could not read meta file '" << metaPath << "'\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        grid_size = static_cast<int>(Nmeta);
        }
    }

    // Alle Ranks erhalten die Metadaten originalSize und grid_size
    MPI_Bcast(&originalSize, 1, MPI_UINT64_MATCHED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&grid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Vorsichtsmaßnahme für zu kleine grid_size (noch kein Fallback implementiert!!!)
    if (nprocs > grid_size) {
        if (rank == 0) std::cerr << "nprocs > grid_size: mindestens ein Rank bekäme 0 Zeilen – für Torus ungültig.\n";
        MPI_Abort(MPI_COMM_WORLD, 2);
    }


    // Lege lokale Zeilen und globalen Offset für jeden Rank fest
    int rows_per_rank = grid_size / nprocs;
    int remainder = grid_size % nprocs;
    int local_rows = rows_per_rank + (rank < remainder ? 1 : 0);
    int offset_rows = rank * rows_per_rank + std::min(rank, remainder);

    // Wird an die Funktionen weitergegeben, um den lokalen Block zu beschreiben
    RowDist dist{};
    dist.grid_size   = grid_size;
    dist.local_rows  = local_rows;
    dist.offset_rows = offset_rows;
    const std::size_t paddedBytes = static_cast<std::size_t>(grid_size) * grid_size;

    // lokaler Speicher für den lokalen Block im Arbeitsspeicher des MPI-Ranks (ohne Halos)
    std::vector<uint8_t> local_core(local_rows * grid_size);

    if (doEncrypt) {
        parallel_read_plain_chunk(inPlain, dist, originalSize, local_core, atomicIO, MPI_COMM_WORLD);
    }
    else {
        parallel_read_cipher_chunk(encBin, dist, paddedBytes, local_core, atomicIO, MPI_COMM_WORLD);
    }
    
    // Flache doppel-Buffer für das Grid (mit Platz für Halo-Zellen)
    std::vector<uint8_t> gridBufA(static_cast<std::size_t>(local_rows + 2) * grid_size);
    std::vector<uint8_t> gridBufB(static_cast<std::size_t>(local_rows + 2) * grid_size);

    // Index-Helfer (i: 0..local_rows+1 inklusive Halos)
    auto idx = [gs = grid_size](int i, int j) -> std::size_t {
        return static_cast<std::size_t>(i) * static_cast<std::size_t>(gs) + static_cast<std::size_t>(j);
    };

    // Kopiere den lokalen Block in das Grid (ohne Halos)
    for (int i = 0; i < local_rows; ++i) {
        std::memcpy(
            gridBufA.data() + idx(i + 1, 0),
            local_core.data() + static_cast<std::size_t>(i) * grid_size,
            static_cast<std::size_t>(grid_size)
        );
    }

    // Rank 0 erzeugt die Wall-Mask oder lädt sie aus der Datei
    Mask wall_mask;
    if (rank == 0) {
        if (doEncrypt) {
            wall_mask = generateRandomWallMask(grid_size, wallDensity, seed);
            saveWallMaskBinary(wall_mask, keyPath.c_str());
            std::cout << "Wall mask generated and saved to " << keyPath << "\n";
        } 
        else {
            wall_mask = loadWallMaskBinary(grid_size, keyPath.c_str());
            std::cout << "Wall mask loaded from " << keyPath << "\n";
        }
    }
    // Alle Ranks erhalten die Wall-Mask
    broadcastMask(wall_mask, MPI_COMM_WORLD);


    // MPI-Topologie für Nachbarschaftskommunikation
    int dims[1]     = { nprocs };   // 1D-Topologie mit nprocs Prozessen
    int periods[1]  = { 1 };        // 1 = periodisch (Torus), 0 = nicht periodisch
    int reorder     = cartReorder;  // Darf MPI die Ranks neu anordnen für bessere Nachbarschaft? Bei Torus nicht erlaubt!!
    MPI_Comm cart_comm; 
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &cart_comm);

    // Rank IDs der Nachbarn (oben, unten) für Halo-Zellen Austausch
    int up, down;
    MPI_Cart_shift(cart_comm, 0, 1, &up, &down);

    // MPI- Requests für Halo-Zellen (4 pro Iteration und pro Grid (2x Receive, 2x Send, 2x Grid))
    MPI_Request halo_reqs_A[4], halo_reqs_B[4];
    // Recvs/Sends für gridA
    MPI_Recv_init(gridBufA.data() + idx(0,             0), grid_size, MPI_BYTE, up,   TAG_NS, cart_comm, &halo_reqs_A[0]);
    MPI_Recv_init(gridBufA.data() + idx(local_rows + 1,0), grid_size, MPI_BYTE, down, TAG_NS, cart_comm, &halo_reqs_A[1]);
    MPI_Send_init(gridBufA.data() + idx(1,             0), grid_size, MPI_BYTE, up,   TAG_NS, cart_comm, &halo_reqs_A[2]);
    MPI_Send_init(gridBufA.data() + idx(local_rows,    0), grid_size, MPI_BYTE, down, TAG_NS, cart_comm, &halo_reqs_A[3]);

    // Recvs/Sends für gridB
    MPI_Recv_init(gridBufB.data() + idx(0,             0), grid_size, MPI_BYTE, up,   TAG_NS, cart_comm, &halo_reqs_B[0]);
    MPI_Recv_init(gridBufB.data() + idx(local_rows + 1,0), grid_size, MPI_BYTE, down, TAG_NS, cart_comm, &halo_reqs_B[1]);
    MPI_Send_init(gridBufB.data() + idx(1,             0), grid_size, MPI_BYTE, up,   TAG_NS, cart_comm, &halo_reqs_B[2]);
    MPI_Send_init(gridBufB.data() + idx(local_rows,    0), grid_size, MPI_BYTE, down, TAG_NS, cart_comm, &halo_reqs_B[3]);

    // Aktive/Passive Buffer- und Request-Sets
    uint8_t*     active_ptr   = gridBufA.data();       
    uint8_t*     target_ptr   = gridBufB.data();         
    MPI_Request* active_reqs  = halo_reqs_A;  
    MPI_Request* idle_reqs    = halo_reqs_B;


    /*------------------------------------------------------------------------------------------------------------------------------------------------------------*/
    // Haupt-Loop: Iterationen der HPP-Regeln
    // 1. Halo-Zellen tauschen (MPI_Send/MPI_Recv)
    // 2. Innenzellen berechnen (mit applyRules)
    // 3. Randzellen berechnen (benötigen Halo-Zellen)
    // 4. Puffer tauschen (grid_current <-> grid_next)
    // 5. Frames speichern (optional)
    /*------------------------------------------------------------------------------------------------------------------------------------------------------------*/

    for (int iter = 0; iter < numIterations; ++iter) {
        // 1) Halo-Transfers des aktuellen Puffers starten
        MPI_Startall(4, active_reqs);

        // 2) Innenbereich berechnen (Zeilen 2..local_rows-1), braucht keine Halos
        if (local_rows >= 3) {
            const int N = grid_size;

            if (doEncrypt) {
                #pragma omp parallel for schedule(static)
                for (int i = 2; i <= local_rows - 1; ++i) {
                    const int gr   = (offset_rows + (i - 1)) % N;       // globale Zeile
                    const int gr_u = (gr - 1 + N) % N;                  // gr-1 (mod N)
                    const int gr_d = (gr + 1) % N;                      // gr+1 (mod N)

                    const uint8_t* __restrict wrow   = wall_mask[gr  ].data();
                    const uint8_t* __restrict wrow_u = wall_mask[gr_u].data();
                    const uint8_t* __restrict wrow_d = wall_mask[gr_d].data();

                    uint8_t* __restrict tgt_row = target_ptr + idx(i, 0);

                    #pragma omp simd
                    for (int j = 0; j < N; ++j) {
                        tgt_row[j] = applyRules_fast<true>(
                            active_ptr, N, i, j, wrow, wrow_u, wrow_d
                        );
                    }
                }
            } else {
                #pragma omp parallel for schedule(static)
                for (int i = 2; i <= local_rows - 1; ++i) {
                    const int gr   = (offset_rows + (i - 1)) % N;
                    const int gr_u = (gr - 1 + N) % N;
                    const int gr_d = (gr + 1) % N;

                    const uint8_t* __restrict wrow   = wall_mask[gr  ].data();
                    const uint8_t* __restrict wrow_u = wall_mask[gr_u].data();
                    const uint8_t* __restrict wrow_d = wall_mask[gr_d].data();

                    uint8_t* __restrict tgt_row = target_ptr + idx(i, 0);

                    #pragma omp simd
                    for (int j = 0; j < N; ++j) {
                        tgt_row[j] = applyRules_fast<false>(
                            active_ptr, N, i, j, wrow, wrow_u, wrow_d
                        );
                    }
                }
            }
        }

        // 3) Auf Halo-Transfers warten, dann Randzeilen (1 und local_rows)
        MPI_Waitall(4, active_reqs, MPI_STATUSES_IGNORE);

        if (local_rows >= 1) {
            const int N    = grid_size;
            const int gr   = (offset_rows + 0) % N;       // i==1 → globale Zeile offset_rows
            const int gr_u = (gr - 1 + N) % N;
            const int gr_d = (gr + 1) % N;

            const uint8_t* __restrict wrow   = wall_mask[gr  ].data();
            const uint8_t* __restrict wrow_u = wall_mask[gr_u].data();
            const uint8_t* __restrict wrow_d = wall_mask[gr_d].data();

            uint8_t* __restrict tgt1 = target_ptr + idx(1, 0);

            if (doEncrypt) {
                #pragma omp simd
                for (int j = 0; j < N; ++j)
                    tgt1[j] = applyRules_fast<true>(active_ptr, N, 1, j, wrow, wrow_u, wrow_d);
            } else {
                #pragma omp simd
                for (int j = 0; j < N; ++j)
                    tgt1[j] = applyRules_fast<false>(active_ptr, N, 1, j, wrow, wrow_u, wrow_d);
            }
        }

        if (local_rows >= 2) {
            const int N    = grid_size;
            const int iL   = local_rows;
            const int gr   = (offset_rows + (iL - 1)) % N;
            const int gr_u = (gr - 1 + N) % N;
            const int gr_d = (gr + 1) % N;

            const uint8_t* __restrict wrow   = wall_mask[gr  ].data();
            const uint8_t* __restrict wrow_u = wall_mask[gr_u].data();
            const uint8_t* __restrict wrow_d = wall_mask[gr_d].data();

            uint8_t* __restrict tgtL = target_ptr + idx(iL, 0);

            if (doEncrypt) {
                #pragma omp simd
                for (int j = 0; j < N; ++j)
                    tgtL[j] = applyRules_fast<true>(active_ptr, N, iL, j, wrow, wrow_u, wrow_d);
            } else {
                #pragma omp simd
                for (int j = 0; j < N; ++j)
                    tgtL[j] = applyRules_fast<false>(active_ptr, N, iL, j, wrow, wrow_u, wrow_d);
            }
        }

        // 4) Puffer und Requests tauschen
        std::swap(active_ptr, target_ptr);
        std::swap(active_reqs, idle_reqs);
    }
    /*
    // Frames speichern
        if (dumpFrames && (iter % frameInterval == 0 || iter == numIterations - 1)) {
            char fname[128];
            std::snprintf(fname, sizeof(fname), "results/frames/frame_%05d.bin", iter);

            std::vector<uint8_t> frame_core(local_rows * grid_size);
            for (int i = 0; i < local_rows; ++i) {
                copy_n_bytes((*active_grid)[i + 1].data(),
                            static_cast<std::size_t>(grid_size),
                            frame_core.data() + i * grid_size);
            }
            dump_frame_parallel(fname, dist, frame_core, paddedBytes, false, MPI_COMM_WORLD);
        }*/

    /*------------------------------------------------------------------------------------------------------------------------------------------------------------*/
    /* Ende des Haupt-Loops */
    /*------------------------------------------------------------------------------------------------------------------------------------------------------------*/


    // Alle Requests und den Topologie-Kommunikator freigeben
    for (auto &r : halo_reqs_A) MPI_Request_free(&r);
    for (auto &r : halo_reqs_B) MPI_Request_free(&r);
    MPI_Comm_free(&cart_comm);


    // Kopiere das Ergebnis (ohne Halo-Zeilen) aus dem aktuellen Puffer in result_core
    std::vector<uint8_t> result_core(static_cast<std::size_t>(local_rows) * grid_size);
    for (int i = 0; i < local_rows; ++i) {
        std::memcpy(
            result_core.data() + static_cast<std::size_t>(i) * grid_size,
            active_ptr        + idx(i + 1, 0),
            static_cast<std::size_t>(grid_size)
        );
    }

    // Ergebnis in Ausgabedatei schreiben (.bin für Encrypt, .txt für Decrypt)
    if (doEncrypt) {
        parallel_write_cipher_chunk(encBin, dist, result_core, paddedBytes, atomicIO, MPI_COMM_WORLD);
        if (rank == 0) {
            write_meta_rank0(metaPath, originalSize, static_cast<uint32_t>(grid_size));
        }
    } 
    else {
        parallel_write_plain_trimmed(outPlain, dist, result_core, originalSize, atomicIO, MPI_COMM_WORLD);
    }
    
    // Zeiterfassung und Ausgabe
    double t_end = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << std::fixed << std::setprecision(6)
                  << "Total runtime: " << (t_end - t_start) << " seconds\n";
    }

    // MPI beenden
    MPI_Finalize();
    return 0;
}

