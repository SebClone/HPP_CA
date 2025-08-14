#pragma once

#include <cstdint>
#include <string>
#include <chrono>

// Betriebsmodus
enum class AppMode : std::uint8_t { Decrypt = 0, Encrypt = 1 };

// Alle Konfigurationen für main.cpp
struct AppConfig {
    AppMode mode = AppMode::Encrypt;
    int iterations = 1000;
    int grid_size = 0;          // 0 => automatisch aus Inputgröße
    int frame_interval = 10;

    std::string input  = "data/message.txt";
    std::string enc_bin = "data/encrypted_full.bin";
    std::string meta    = "data/encrypted_full.meta";
    std::string key     = "data/wall_mask.key";
    std::string output  = "data/decrypted_message.txt";

    bool dump_frames = false;    // true = Frames speichern
    double wall_density = 0.10;  // nur im Encrypt-Fall relevant
    std::uint64_t seed = 0;      // 0 => nicht-deterministisch (zeitbasiert)
};

namespace appcfg {
    inline constexpr AppMode mode = AppMode::Encrypt;
    inline constexpr int iterations = 1000;
    inline constexpr int grid_size = 0;
    inline constexpr int frame_interval = 10;

    inline constexpr const char* input   = "data/message.txt";
    inline constexpr const char* enc_bin = "data/encrypted_full.bin";
    inline constexpr const char* meta    = "data/encrypted_full.meta";
    inline constexpr const char* key     = "data/wall_mask.key";
    inline constexpr const char* output  = "data/decrypted_message.txt";

    inline constexpr bool dump_frames  = false;

    inline constexpr double wall_density = 0.10;
    inline constexpr std::uint64_t seed = 0ULL;
} 


inline AppConfig get_config() {
    AppConfig c;
    c.mode           = appcfg::mode;
    c.iterations     = appcfg::iterations;
    c.grid_size      = appcfg::grid_size;
    c.frame_interval = appcfg::frame_interval;

    c.input   = appcfg::input;
    c.enc_bin = appcfg::enc_bin;
    c.meta    = appcfg::meta;
    c.key     = appcfg::key;
    c.output  = appcfg::output;

    c.dump_frames  = appcfg::dump_frames;

    c.wall_density = appcfg::wall_density;
    c.seed         = appcfg::seed;

    if (c.seed == 0ULL) {
        const auto now = static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
        c.seed = now;
    }
    return c;
}

inline bool validate_config(const AppConfig& cfg, std::string& error_msg) {
    if (cfg.iterations <= 0) { error_msg = "iterations must be > 0"; return false; }
    if (cfg.grid_size < 0)   { error_msg = "grid-size must be >= 0"; return false; }
    if (cfg.wall_density < 0.0 || cfg.wall_density > 1.0) {
        error_msg = "wall-density must be in [0,1]"; return false;
    }
    if (cfg.frame_interval <= 0) { error_msg = "frame-interval must be > 0"; return false; }
    (void)cfg;
    return true;
}