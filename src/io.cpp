#include "io.hpp"   

#include <algorithm>
#include <fstream>
#include <iostream>

namespace {

// kleine Helfer
inline MPI_Offset byte_offset(const RowDist& d) {
    return static_cast<MPI_Offset>(d.offset_rows) * d.grid_size; // 1 Byte/Elem
}
inline std::size_t local_bytes(const RowDist& d) {
    return static_cast<std::size_t>(d.local_rows) * d.grid_size;
}

// einheitliche MPI-Fehlerprüfung
inline void mpi_check(int rc, const char* where, MPI_Comm comm) {
    if (rc == MPI_SUCCESS) return;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    char errstr[MPI_MAX_ERROR_STRING]; int len=0;
    MPI_Error_string(rc, errstr, &len);
    std::cerr << "[r" << rank << "] MPI error at " << where << ": "
              << std::string(errstr, len) << "\n";
    MPI_Abort(comm, 1);
}

} // namespace

// ---------- Lesen/Schreiben: Plain -----------

void parallel_read_plain_chunk(const std::string& path,
                               const RowDist& dist,
                               std::uint64_t global_plain_size,
                               std::vector<std::uint8_t>& out,
                               bool atomic_io,
                               MPI_Comm comm) {
    out.assign(local_bytes(dist), 0); // init mit 0 (Padding)

    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh), "open plain", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity plain", comm);

    const MPI_Offset base = byte_offset(dist);
    const std::size_t want = local_bytes(dist);

    std::size_t to_read = 0;
    if (static_cast<std::uint64_t>(base) < global_plain_size) {
        const auto remain = global_plain_size - static_cast<std::uint64_t>(base);
        to_read = static_cast<std::size_t>(std::min<std::uint64_t>(want, remain));
    }

    if (to_read > 0) {
        mpi_check(MPI_File_read_at_all(fh, base, out.data(),
                                       static_cast<int>(to_read), MPI_BYTE, MPI_STATUS_IGNORE),
                  "read_at_all plain", comm);
    }
    mpi_check(MPI_File_close(&fh), "close plain", comm);
}

void parallel_read_cipher_chunk(const std::string& path,
                                const RowDist& dist,
                                std::size_t padded_size_bytes,
                                std::vector<std::uint8_t>& out,
                                bool atomic_io,
                                MPI_Comm comm) {
    out.resize(local_bytes(dist));

    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh), "open cipher", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity cipher", comm);

    MPI_Offset fsz = 0;
    mpi_check(MPI_File_get_size(fh, &fsz), "get_size cipher", comm);
    if (static_cast<std::size_t>(fsz) != padded_size_bytes) {
        int rank=0; MPI_Comm_rank(comm, &rank);
        if (rank==0) {
            std::cerr << "ERROR: '" << path << "' has " << fsz
                      << " bytes, expected " << padded_size_bytes << "\n";
        }
        MPI_File_close(&fh);
        MPI_Abort(comm, 1);
    }

    mpi_check(MPI_File_read_at_all(fh, byte_offset(dist), out.data(),
                                   static_cast<int>(out.size()), MPI_BYTE, MPI_STATUS_IGNORE),
              "read_at_all cipher", comm);
    mpi_check(MPI_File_close(&fh), "close cipher", comm);
}

void parallel_write_cipher_chunk(const std::string& path,
                                 const RowDist& dist,
                                 const std::vector<std::uint8_t>& data,
                                 std::size_t padded_size_bytes,
                                 bool atomic_io,
                                 MPI_Comm comm) {
    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(),
                            MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &fh),
              "open cipher write", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity cipher write", comm);

    int rank=0; MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        mpi_check(MPI_File_set_size(fh, static_cast<MPI_Offset>(padded_size_bytes)),
                  "set_size cipher write", comm);
    }
    MPI_Barrier(comm);

    mpi_check(MPI_File_write_at_all(fh, byte_offset(dist), data.data(),
                                    static_cast<int>(data.size()),
                                    MPI_BYTE, MPI_STATUS_IGNORE),
              "write_at_all cipher", comm);
    mpi_check(MPI_File_close(&fh), "close cipher write", comm);
}

void parallel_write_plain_trimmed(const std::string& path,
                                  const RowDist& dist,
                                  const std::vector<std::uint8_t>& data,
                                  std::uint64_t original_size,
                                  bool atomic_io,
                                  MPI_Comm comm) {
    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(),
                            MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &fh),
              "open plain write", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity plain write", comm);

    int rank=0; MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        mpi_check(MPI_File_set_size(fh, static_cast<MPI_Offset>(original_size)),
                  "set_size plain write", comm);
    }
    MPI_Barrier(comm);

    const MPI_Offset base = byte_offset(dist);
    const std::size_t have = data.size();

    std::size_t to_write = 0;
    if (static_cast<std::uint64_t>(base) < original_size) {
        const auto remain = original_size - static_cast<std::uint64_t>(base);
        to_write = static_cast<std::size_t>(std::min<std::uint64_t>(have, remain));
    }

    if (to_write > 0) {
        mpi_check(MPI_File_write_at_all(fh, base, data.data(),
                                        static_cast<int>(to_write),
                                        MPI_BYTE, MPI_STATUS_IGNORE),
                  "write_at_all plain", comm);
    }
    mpi_check(MPI_File_close(&fh), "close plain write", comm);
}

void dump_frame_parallel(const std::string& path,
                         const RowDist& dist,
                         const std::vector<std::uint8_t>& data,
                         std::size_t padded_size_bytes,
                         bool atomic_io,
                         MPI_Comm comm) {
    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(),
                            MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &fh),
              "open frame", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity frame", comm);

    int rank=0; MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        mpi_check(MPI_File_set_size(fh, static_cast<MPI_Offset>(padded_size_bytes)),
                  "set_size frame", comm);
    }
    MPI_Barrier(comm);

    mpi_check(MPI_File_write_at_all(fh, byte_offset(dist), data.data(),
                                    static_cast<int>(data.size()),
                                    MPI_BYTE, MPI_STATUS_IGNORE),
              "write_at_all frame", comm);
    mpi_check(MPI_File_close(&fh), "close frame", comm);
}

// ---------- Meta (kleine Files, Rank 0) ----------

bool read_meta_rank0(const std::string& meta_path,
                     std::uint64_t& original_size,
                     std::uint32_t& grid_size_n) {
    std::ifstream in(meta_path, std::ios::binary);
    if (!in) return false;
    in.read(reinterpret_cast<char*>(&original_size), sizeof(original_size));
    in.read(reinterpret_cast<char*>(&grid_size_n),   sizeof(grid_size_n));
    return static_cast<bool>(in);
}

bool write_meta_rank0(const std::string& meta_path,
                      std::uint64_t original_size,
                      std::uint32_t grid_size_n) {
    std::ofstream out(meta_path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out.write(reinterpret_cast<const char*>(&original_size), sizeof(original_size));
    out.write(reinterpret_cast<const char*>(&grid_size_n),   sizeof(grid_size_n));
    return static_cast<bool>(out);
}