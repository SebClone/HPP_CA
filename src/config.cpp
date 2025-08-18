#include "config.hpp"
#include <chrono>
#include <string>

AppConfig get_config()
{
    AppConfig c;
    c.mode = appcfg::mode;
    c.iterations = appcfg::iterations;
    c.grid_size = appcfg::grid_size;
    c.frame_interval = appcfg::frame_interval;

    // Festlegen, ob Bild oder Text verarbeitet wird
    // Text
    c.input = appcfg::input_text;
    c.output = appcfg::output_text;

    // // Bild
    // c.input = appcfg::input_img;
    // c.output = appcfg::output_img;

    c.enc_bin = appcfg::enc_bin;
    c.meta = appcfg::meta;
    c.key = appcfg::key;

    c.dump_frames = appcfg::dump_frames;
    c.wall_density = appcfg::wall_density;
    c.seed = appcfg::seed;

    if (c.seed == 0ULL)
    {
        const auto now = static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
        c.seed = now;
    }
    return c;
}

/**
 * @brief Validates the application configuration.
 *
 * Checks if the configuration parameters are within valid ranges.
 * @param cfg The application configuration to validate.
 * @param error_msg Output parameter to store error messages if validation fails.
 * @return true if the configuration is valid, false otherwise.
 */
bool validate_config(const AppConfig &cfg, std::string &error_msg)
{
    if (cfg.iterations <= 0)
    {
        error_msg = "iterations must be > 0";
        return false;
    }
    if (cfg.grid_size < 0)
    {
        error_msg = "grid-size must be >= 0";
        return false;
    }
    if (cfg.wall_density < 0.0 || cfg.wall_density > 1.0)
    {
        error_msg = "wall-density must be in [0,1]";
        return false;
    }
    if (cfg.frame_interval <= 0)
    {
        error_msg = "frame-interval must be > 0";
        return false;
    }
    (void)cfg;
    return true;
}

/**
 * @brief Parses the application mode from command line arguments.
 *
 * Checks for the presence of specific flags to determine the mode (Encrypt or Decrypt).
 * For example, "--mode encrypt" or "--encrypt" for encryption mode, and "--mode decrypt" or "--decrypt" for decryption mode.
 * If no flags are provided, it returns the default mode.
 * @param argc The number of command line arguments.
 * @param argv The array of command line arguments.
 * @param deflt The default mode if no flags are provided.
 * @return The parsed application mode.
 */
AppMode parse_mode_from_cli(int argc, char **argv, AppMode deflt)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--mode" || a == "-m")
        {
            if (i + 1 < argc)
            {
                const std::string v = argv[++i];
                if (v == "encrypt" || v == "enc" || v == "e")
                    return AppMode::Encrypt;
                if (v == "decrypt" || v == "dec" || v == "d")
                    return AppMode::Decrypt;
            }
        }
        else if (a == "--encrypt" || a == "-E")
        {
            return AppMode::Encrypt;
        }
        else if (a == "--decrypt" || a == "-D")
        {
            return AppMode::Decrypt;
        }
    }
    return deflt;
}

/**
 * @brief Parses the grid size from command line arguments.
 *
 * Checks for the presence of the "--grid" or "-g" flag to set the grid size.
 * If the flag is not present, it returns the fallback value.
 * @param argc The number of command line arguments.
 * @param argv The array of command line arguments.
 * @param fallback The fallback value for grid size if no flag is provided.
 * @return The parsed grid size.
 */
int parse_grid_from_cli(int argc, char **argv, int fallback)
{
    int N = fallback;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if ((a == "--grid" || a == "-g") && i + 1 < argc)
        {
            try
            {
                N = std::stoi(argv[++i]);
            }
            catch (...)
            {
                // ignore
            }
        }
    }
    return N;
}

/**
 * @brief Parses the number of iterations from command line arguments.
 *
 * Checks for the presence of the "--iters" or "-I" flag to set the number of iterations.
 * If the flag is not present, it returns the fallback value.
 * @param argc The number of command line arguments.
 * @param argv The array of command line arguments.
 * @param fallback The fallback value for iterations if no flag is provided.
 * @return The parsed number of iterations.
 */
int parse_iters_from_cli(int argc, char **argv, int fallback)
{
    int iters = fallback;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if ((a == "--iters" || a == "-I") && i + 1 < argc)
        {
            try
            {
                iters = std::stoi(argv[++i]);
            }
            catch (...)
            {
                // ignore
            }
        }
    }
    return iters;
}