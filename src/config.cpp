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

    c.input = appcfg::input;
    c.enc_bin = appcfg::enc_bin;
    c.meta = appcfg::meta;
    c.key = appcfg::key;
    c.output = appcfg::output;

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
                // ignore, keep fallback
            }
        }
    }
    return N;
}

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