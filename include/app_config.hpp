#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <mpi.h>

enum class AppMode : uint8_t { Decrypt = 0, Encrypt = 1 };

struct AppConfig {
    AppMode mode = AppMode::Decrypt;
    int iterations = 1000;
    int grid_size = 0;          // 0 => automatisch aus Inputgröße
    int frame_interval = 5;
    int verbosity = 1;          // 0=quiet,1=normal,2=verbose

    std::string input  = "message.txt";
    std::string enc_bin = "encrypted_full.bin";
    std::string meta    = "encrypted_full.meta";
    std::string key     = "wall_mask.key";
    std::string output  = "decrypted_message.txt";  

    bool dump_frames = true;
    bool atomic_io = false;      // MPI_File_set_atomicity
    bool cart_reorder = false;   // MPI_Cart_create reorder

    double wall_density = 0.10;  // nur im Encrypt-Fall relevant
    std::uint64_t seed = 0;      // 0 => non-deterministisch (z.B. Zeitstempel)
};

// ---- API ----
AppConfig parse_args_or_default(int argc, char** argv, int mpi_rank);

bool validate_config(const AppConfig& cfg, std::string& error_msg);

void bcast_config(AppConfig& cfg, int root, MPI_Comm comm);

void print_usage(const char* progname);

inline bool is_verbose(const AppConfig& cfg) { return cfg.verbosity >= 2; }
inline bool is_quiet(const AppConfig& cfg)   { return cfg.verbosity == 0; }