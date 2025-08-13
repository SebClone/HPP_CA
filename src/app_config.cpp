#include "app_config.hpp"      

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <type_traits>

// ---------- CLI ----------
static void assign_if_set_int(const char* optarg_c, int& out) {
    if (!optarg_c) return;
    char* endp = nullptr;
    long v = std::strtol(optarg_c, &endp, 10);
    if (endp && *endp == '\0') out = static_cast<int>(v);
}

static void assign_if_set_u64(const char* optarg_c, std::uint64_t& out) {
    if (!optarg_c) return;
    char* endp = nullptr;
    unsigned long long v = std::strtoull(optarg_c, &endp, 10);
    if (endp && *endp == '\0') out = static_cast<std::uint64_t>(v);
}

static void assign_if_set_double(const char* optarg_c, double& out) {
    if (!optarg_c) return;
    char* endp = nullptr;
    double v = std::strtod(optarg_c, &endp);
    if (endp && *endp == '\0') out = v;
}

void print_usage(const char* progname) {
    std::fprintf(stderr,
R"(Usage: %s [--encrypt | --decrypt]
       [--input file] [--enc-bin file] [--meta file] [--key file] [--output file]
       [--iterations N] [--grid-size N] [--wall-density f] [--seed S]
       [--dump-frames 0|1] [--frame-interval K]
       [--atomic-io 0|1] [--reorder 0|1]
       [-v | -q]

Modes (choose one; default: --decrypt):
  --encrypt                 Verschlüsseln
  --decrypt                 Entschlüsseln

I/O:
  --input PATH              Klartext-Eingabe (Encrypt) [default: message.txt]
  --enc-bin PATH            Cipher-Binärdatei [default: encrypted_full.bin]
  --meta PATH               Meta-Datei [default: encrypted_full.meta]
  --key PATH                Wandmaske [default: wall_mask.key]
  --output PATH             Entschlüsselter Output [default: decrypted_message.txt]

Rechnen:
  --iterations N            Iterationen (default 1000)
  --grid-size N             optionales Override der Rastergröße
  --wall-density f          Dichte der Wände (Encrypt) (default 0.10)
  --seed S                  RNG Seed (0 => auto)

Debug/Perf:
  --dump-frames 0|1         Frames schreiben (default 1)
  --frame-interval K        Intervall für Frames (default 5)
  --atomic-io 0|1           MPI atomic I/O (default 0)
  --reorder 0|1             MPI_Cart_create reorder (default 0)
  -v                        verbose (mehr Logs auf Rank 0)
  -q                        quiet

Examples:
  %s --encrypt --input msg.txt --enc-bin out.bin --meta out.meta --key k.bin
  %s --decrypt --enc-bin out.bin --meta out.meta --key k.bin --output clear.txt
)",
    progname, progname, progname);
}

AppConfig parse_args_or_default(int argc, char** argv, int mpi_rank) {
    AppConfig cfg;

    if (mpi_rank != 0) {
        // Non-Root: nichts parsen; Root broadcastet später den echten cfg.
        return cfg;
    }

    static option long_opts[] = {
        {"encrypt",       no_argument,       nullptr,  1},
        {"decrypt",       no_argument,       nullptr,  2},
        {"input",         required_argument, nullptr,  3},
        {"enc-bin",       required_argument, nullptr,  4},
        {"meta",          required_argument, nullptr,  5},
        {"key",           required_argument, nullptr,  6},
        {"output",        required_argument, nullptr,  7},
        {"iterations",    required_argument, nullptr,  8},
        {"grid-size",     required_argument, nullptr,  9},
        {"wall-density",  required_argument, nullptr, 10},
        {"seed",          required_argument, nullptr, 11},
        {"dump-frames",   required_argument, nullptr, 12},
        {"frame-interval",required_argument, nullptr, 13},
        {"atomic-io",     required_argument, nullptr, 14},
        {"reorder",       required_argument, nullptr, 15},
        {nullptr,         0,                 nullptr,  0}
    };

    int opt, long_idx;
    while ((opt = getopt_long(argc, argv, "vq", long_opts, &long_idx)) != -1) {
        switch (opt) {
            case 0: break; // sollte nicht auftreten
            case 1: cfg.mode = AppMode::Encrypt; break;
            case 2: cfg.mode = AppMode::Decrypt; break;
            case 3: if (optarg) cfg.input = optarg; break;
            case 4: if (optarg) cfg.enc_bin = optarg; break;
            case 5: if (optarg) cfg.meta    = optarg; break;
            case 6: if (optarg) cfg.key     = optarg; break;
            case 7: if (optarg) cfg.output  = optarg; break;
            case 8: assign_if_set_int(optarg, cfg.iterations); break;
            case 9: assign_if_set_int(optarg, cfg.grid_size); break;
            case 10: assign_if_set_double(optarg, cfg.wall_density); break;
            case 11: assign_if_set_u64(optarg, cfg.seed); break;
            case 12: { int v=cfg.dump_frames; assign_if_set_int(optarg, v); cfg.dump_frames = (v!=0); } break;
            case 13: assign_if_set_int(optarg, cfg.frame_interval); break;
            case 14: { int v=cfg.atomic_io; assign_if_set_int(optarg, v); cfg.atomic_io = (v!=0); } break;
            case 15: { int v=cfg.cart_reorder; assign_if_set_int(optarg, v); cfg.cart_reorder = (v!=0); } break;
            case 'v': cfg.verbosity = 2; break;
            case 'q': cfg.verbosity = 0; break;
            default:
                print_usage(argv[0]);
                std::exit(2);
        }
    }

    if (cfg.seed == 0) {
        auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        cfg.seed = static_cast<std::uint64_t>(now);
    }

    // Minimaler Plausibilitätscheck auf Rank 0; harte Fehler lassen wir den Caller handeln.
    std::string err;
    if (!validate_config(cfg, err)) {
        std::cerr << "[r0] config error: " << err << "\n\n";
        print_usage(argv[0]);
        std::exit(2);
    }
    return cfg;
}

bool validate_config(const AppConfig& cfg, std::string& error_msg) {
    if (cfg.iterations <= 0) {
        error_msg = "iterations must be > 0";
        return false;
    }
    if (cfg.grid_size < 0) {
        error_msg = "grid-size must be >= 0";
        return false;
    }
    if (cfg.wall_density < 0.0 || cfg.wall_density > 1.0) {
        error_msg = "wall-density must be in [0,1]";
        return false;
    }
    if (cfg.frame_interval <= 0) {
        error_msg = "frame-interval must be > 0";
        return false;
    }
    // Dateien: keine harte Prüfung hier; die I/O-Schicht prüft existence/size.
    (void)cfg;
    return true;
}

// ---------- Serialization für MPI_Bcast ----------

namespace {

template <class T>
void push_pod(std::vector<std::uint8_t>& buf, const T& v) {
    static_assert(std::is_trivially_copyable<T>::value, "POD only");
    const std::uint8_t* p = reinterpret_cast<const std::uint8_t*>(&v);
    buf.insert(buf.end(), p, p + sizeof(T));
}

template <class T>
T read_pod(const std::vector<std::uint8_t>& buf, size_t& off) {
    T v{};
    std::memcpy(&v, buf.data() + off, sizeof(T));
    off += sizeof(T);
    return v;
}

void push_string(std::vector<std::uint8_t>& buf, const std::string& s) {
    std::uint32_t n = static_cast<std::uint32_t>(s.size());
    push_pod(buf, n);
    buf.insert(buf.end(), s.begin(), s.end());
}

std::string read_string(const std::vector<std::uint8_t>& buf, size_t& off) {
    std::uint32_t n = read_pod<std::uint32_t>(buf, off);
    std::string s;
    s.resize(n);
    if (n) std::memcpy(&s[0], buf.data() + off, n);
    off += n;
    return s;
}

} // namespace

static std::vector<std::uint8_t> serialize_cfg(const AppConfig& c) {
    std::vector<std::uint8_t> b;
    b.reserve(256);

    // Ordnung beibehalten!
    push_pod(b, c.mode);
    push_pod(b, c.iterations);
    push_pod(b, c.grid_size);
    push_pod(b, c.frame_interval);
    push_pod(b, c.verbosity);

    push_string(b, c.input);
    push_string(b, c.enc_bin);
    push_string(b, c.meta);
    push_string(b, c.key);
    push_string(b, c.output);

    push_pod(b, c.dump_frames);
    push_pod(b, c.atomic_io);
    push_pod(b, c.cart_reorder);

    push_pod(b, c.wall_density);
    push_pod(b, c.seed);

    return b;
}

static AppConfig deserialize_cfg(const std::vector<std::uint8_t>& b) {
    AppConfig c;
    size_t off = 0;

    c.mode          = read_pod<AppMode>(b, off);
    c.iterations    = read_pod<int>(b, off);
    c.grid_size     = read_pod<int>(b, off);
    c.frame_interval= read_pod<int>(b, off);
    c.verbosity     = read_pod<int>(b, off);

    c.input   = read_string(b, off);
    c.enc_bin = read_string(b, off);
    c.meta    = read_string(b, off);
    c.key     = read_string(b, off);
    c.output  = read_string(b, off);

    c.dump_frames  = read_pod<bool>(b, off);
    c.atomic_io    = read_pod<bool>(b, off);
    c.cart_reorder = read_pod<bool>(b, off);

    c.wall_density = read_pod<double>(b, off);
    c.seed         = read_pod<std::uint64_t>(b, off);

    return c;
}

void bcast_config(AppConfig& cfg, int root, MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    std::vector<std::uint8_t> payload;

    if (rank == root) {
        payload = serialize_cfg(cfg);
    }
    // erst Länge senden
    std::uint32_t n = static_cast<std::uint32_t>(payload.size());
    MPI_Bcast(&n, 1, MPI_UNSIGNED, root, comm);

    if (rank != root) {
        payload.resize(n);
    }

    if (n > 0) {
        MPI_Bcast(payload.data(), static_cast<int>(n), MPI_BYTE, root, comm);
    }

    if (rank != root) {
        cfg = deserialize_cfg(payload);
    }
}
