#pragma once
#include <cstdint>
#include <string>

// Betriebsmodus
enum class AppMode : std::uint8_t
{
    Decrypt = 0,
    Encrypt = 1
};

// Alle Konfigurationen für main.cpp
struct AppConfig
{
    AppMode mode = AppMode::Encrypt;
    int iterations = 1000;
    int grid_size = 0; // 0 => automatisch aus Inputgröße
    int frame_interval = 10;

    std::string input = "data/message.txt";
    std::string enc_bin = "data/encrypted_full.bin";
    std::string meta = "data/encrypted_full.meta";
    std::string key = "data/wall_mask.key";
    std::string output = "data/decrypted_message.txt";

    bool dump_frames = false;   // true = Frames speichern
    double wall_density = 0.10; // nur im Encrypt-Fall relevant
    std::uint64_t seed = 0;     // 0 => nicht-deterministisch (zeitbasiert)
};

// Defaultwerte (können im Header als inline constexpr bleiben)
namespace appcfg
{
    inline constexpr AppMode mode = AppMode::Encrypt;
    inline constexpr int iterations = 1000;
    inline constexpr int grid_size = 0;
    inline constexpr int frame_interval = 10;

    inline constexpr const char *input = "data/message.txt";
    inline constexpr const char *enc_bin = "data/encrypted_full.bin";
    inline constexpr const char *meta = "data/encrypted_full.meta";
    inline constexpr const char *key = "data/wall_mask.key";
    inline constexpr const char *output = "data/decrypted_message.txt";

    inline constexpr bool dump_frames = false;
    inline constexpr double wall_density = 0.10;
    inline constexpr std::uint64_t seed = 0ULL;
}

// Deklerationen der Konfigurationsfunktionen
AppConfig get_config();
bool validate_config(const AppConfig &cfg, std::string &error_msg);

AppMode parse_mode_from_cli(int argc, char **argv, AppMode deflt);
int parse_grid_from_cli(int argc, char **argv, int fallback);
int parse_iters_from_cli(int argc, char **argv, int fallback);