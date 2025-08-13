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
    RowDist dist{ .grid_size = grid_size, .local_rows = local_rows, .offset_rows = offset_rows };
    const std::size_t paddedBytes = static_cast<std::size_t>(grid_size) * grid_size;

    // lokaler Speicher für den lokalen Block im Arbeitsspeicher des MPI-Ranks (ohne Halos)
    std::vector<uint8_t> local_core(local_rows * grid_size);

    if (doEncrypt) {
        parallel_read_plain_chunk(inPlain, dist, originalSize, local_core, atomicIO, MPI_COMM_WORLD);
    }
    else {
        parallel_read_cipher_chunk(encBin, dist, paddedBytes, local_core, /*atomic*/false, MPI_COMM_WORLD);
    }
    
    // doppel-Buffer des Grids mit Platz für Halo-Zellen 
    Matrix grid_current(local_rows + 2, std::vector<uint8_t>(grid_size));
    Matrix grid_next = grid_current;

    // Kopiere lokale Daten in den Grid-Buffer
    for (int i = 0; i < local_rows; ++i) {
        copy_n_bytes(
            local_core.data() + i * grid_size,
            static_cast<std::size_t>(grid_size),
            grid_current[i + 1].data());
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
    MPI_Request halo_reqs_current[4], halo_reqs_next[4];
    MPI_Recv_init(grid_current[0].data(),               grid_size, MPI_BYTE, up,    TAG_NS, cart_comm, &halo_reqs_current[0]);
    MPI_Recv_init(grid_current[local_rows + 1].data(),  grid_size, MPI_BYTE, down,  TAG_NS, cart_comm, &halo_reqs_current[1]);
    MPI_Send_init(grid_current[1].data(),               grid_size, MPI_BYTE, up,    TAG_NS, cart_comm, &halo_reqs_current[2]);
    MPI_Send_init(grid_current[local_rows].data(),      grid_size, MPI_BYTE, down,  TAG_NS, cart_comm, &halo_reqs_current[3]);
    MPI_Recv_init(grid_next[0].data(),                  grid_size, MPI_BYTE, up,    TAG_NS, cart_comm, &halo_reqs_next[0]);
    MPI_Recv_init(grid_next[local_rows + 1].data(),     grid_size, MPI_BYTE, down,  TAG_NS, cart_comm, &halo_reqs_next[1]);
    MPI_Send_init(grid_next[1].data(),                  grid_size, MPI_BYTE, up,    TAG_NS, cart_comm, &halo_reqs_next[2]);
    MPI_Send_init(grid_next[local_rows].data(),         grid_size, MPI_BYTE, down,  TAG_NS, cart_comm, &halo_reqs_next[3]);

    Matrix* active_grid = &grid_current;            // Quelle zum lesen (aktuelles Grid)
    Matrix* target_grid = &grid_next;               // Ziel zum schreiben (nächstes Grid)
    MPI_Request* active_reqs = halo_reqs_current;   // Aktiver Satz persistenter Halo-Requests (gehört zu grid_current / active_grid)
    MPI_Request* inactive_reqs = halo_reqs_next;    // Anderer Satz Requests (für grid_next / target_grid). Wird nach jeder Iteration mit active_reqs geswapped

    for (int iter = 0; iter < numIterations; ++iter) {
        MPI_Startall(4, active_reqs);

        // Innenbereich berechnen (benötigt keine Halo-Zellen)
        if (local_rows >= 3) {
            for (int i = 2; i <= local_rows - 1; i++) {
                for (int j = 0; j < grid_size; ++j) {
                    (*target_grid)[i][j] = applyRules(*active_grid, wall_mask, doEncrypt, i, j, offset_rows);
                }
            }
        }

        MPI_Waitall(4, active_reqs, MPI_STATUSES_IGNORE);

        // Zellen im Randbereich berechnen (benötigen Halo-Zellen)
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
        // Puffer und Requests tauschen
        std::swap(active_grid, target_grid);
        std::swap(active_reqs, inactive_reqs);

        // Frames speichern
        if (dumpFrames && (iter % frameInterval == 0 || iter == numIterations - 1)) {
            char fname[128];
            std::snprintf(fname, sizeof(fname), "frame_%05d.bin", iter);

            std::vector<uint8_t> frame_core(local_rows * grid_size);
            for (int i = 0; i < local_rows; ++i) {
                copy_n_bytes((*active_grid)[i + 1].data(),
                            static_cast<std::size_t>(grid_size),
                            frame_core.data() + i * grid_size);
            }
            dump_frame_parallel(fname, dist, frame_core, paddedBytes, /*atomic*/false, MPI_COMM_WORLD);
        }
    }

    // Alle RRequests und den Topologie-Kommunikator freigeben
    for (auto &r : halo_reqs_current) MPI_Request_free(&r);
    for (auto &r : halo_reqs_next)    MPI_Request_free(&r);
    MPI_Comm_free(&cart_comm);

    // Kopiere das Ergebnis in den lokalen Speicher (ohne Halo-Zellen)
    std::vector<uint8_t> result_core(local_rows * grid_size);
    for (int i = 0; i < local_rows; ++i)
    {
        copy_n_bytes(
            (*active_grid)[i + 1].data(),
            static_cast<std::size_t>(grid_size),
            result_core.data() + i * grid_size);
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

