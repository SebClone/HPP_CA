#include "io.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>

namespace
{
    inline MPI_Offset byte_offset(const RowDist &d)
    {
        return static_cast<MPI_Offset>(d.offset_rows) * d.grid_size; // 1 Byte/Elem
    }
    inline std::size_t local_bytes(const RowDist &d)
    {
        return static_cast<std::size_t>(d.local_rows) * d.grid_size;
    }

    // einheitliche MPI-Fehlerprüfung
    inline void mpi_check(int rc, const char *where, MPI_Comm comm)
    {
        if (rc == MPI_SUCCESS)
            return;
        int rank = 0;
        MPI_Comm_rank(comm, &rank);
        char errstr[MPI_MAX_ERROR_STRING];
        int len = 0;
        MPI_Error_string(rc, errstr, &len);
        std::cerr << "[r" << rank << "] MPI error at " << where << ": "
                  << std::string(errstr, len) << "\n";
        MPI_Abort(comm, 1);
    }
}

void parallel_read_plain_chunk(const std::string &path,
                               const RowDist &dist,
                               std::uint64_t original_size,
                               std::uint64_t start_offset,
                               std::vector<std::uint8_t> &out,
                               bool atomic_io,
                               MPI_Comm comm)
{
    const std::uint64_t base = (std::uint64_t)dist.offset_rows * dist.grid_size;
    const std::uint64_t len = (std::uint64_t)dist.local_rows * dist.grid_size;
    const std::uint64_t end = base + len;

    // Nachrichtenbereich [start_offset, start_offset + original_size)
    const std::uint64_t msgL = start_offset;
    const std::uint64_t msgR = start_offset + original_size;

    // Overlap [L,R)
    const std::uint64_t L = std::max(base, msgL);
    const std::uint64_t R = std::min(end, msgR);

    out.assign((std::size_t)len, 0); // Pre-Fill mit 0 (Padding)

    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh), "open plain", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity plain", comm);

    if (R > L)
    {
        const std::uint64_t read_bytes = R - L;
        const MPI_Offset file_off = (MPI_Offset)(L - msgL);   // Input-File startet bei 0
        const std::size_t dest_off = (std::size_t)(L - base); // Position im lokalen Puffer
        mpi_check(MPI_File_read_at(fh, file_off,
                                   out.data() + dest_off,
                                   (int)read_bytes, MPI_BYTE, MPI_STATUS_IGNORE),
                  "read_at plain(centered)", comm);
    }

    mpi_check(MPI_File_close(&fh), "close plain", comm);
}
// == READ (Cipher) ==
void parallel_read_cipher_chunk(const std::string &path,
                                const RowDist &dist,
                                std::size_t padded_size_bytes,
                                std::vector<std::uint8_t> &out,
                                bool atomic_io,
                                MPI_Comm comm)
{
    out.resize(static_cast<std::size_t>(dist.local_rows) * dist.grid_size);

    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh),
              "open cipher", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0),
              "set_atomicity cipher", comm);

    MPI_Offset fsz = 0;
    mpi_check(MPI_File_get_size(fh, &fsz), "get_size cipher", comm);
    if (static_cast<std::size_t>(fsz) != padded_size_bytes)
    {
        int r = 0;
        MPI_Comm_rank(comm, &r);
        if (r == 0)
        {
            std::cerr << "ERROR: '" << path << "' has " << fsz
                      << " bytes, expected " << padded_size_bytes << "\n";
        }
        MPI_File_close(&fh);
        MPI_Abort(comm, 1);
    }

    // Leichte Synchronisation: alle Ranks kennen die erwartete Dateigröße
    std::size_t expect = padded_size_bytes; // broadcast zur Sicherheit
    MPI_Bcast(&expect, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

    mpi_check(MPI_File_read_at(fh,
                               static_cast<MPI_Offset>(dist.offset_rows) * dist.grid_size,
                               out.data(),
                               static_cast<int>(out.size()),
                               MPI_BYTE, MPI_STATUS_IGNORE),
              "read_at cipher", comm);

    mpi_check(MPI_File_close(&fh), "close cipher", comm);
}
// == WRITE (Cipher) ==
void parallel_write_cipher_chunk(const std::string &path,
                                 const RowDist &dist,
                                 const std::vector<std::uint8_t> &data,
                                 std::size_t padded_size_bytes,
                                 bool atomic_io,
                                 MPI_Comm comm)
{
    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(),
                            MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &fh),
              "open cipher write", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0),
              "set_atomicity cipher write", comm);

    mpi_check(MPI_File_write_at(fh,
                                static_cast<MPI_Offset>(dist.offset_rows) * dist.grid_size,
                                data.data(),
                                static_cast<int>(data.size()),
                                MPI_BYTE, MPI_STATUS_IGNORE),
              "write_at cipher", comm);

    mpi_check(MPI_File_close(&fh), "close cipher write", comm);

    int r0 = 0;
    MPI_Comm_rank(comm, &r0);
    if (r0 == 0)
    {
        MPI_File f2;
        mpi_check(MPI_File_open(MPI_COMM_SELF, path.c_str(),
                                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f2),
                  "open cipher resize", MPI_COMM_SELF);
        mpi_check(MPI_File_set_size(f2, static_cast<MPI_Offset>(padded_size_bytes)),
                  "set_size cipher post", MPI_COMM_SELF);
        mpi_check(MPI_File_close(&f2), "close cipher resize", MPI_COMM_SELF);
    }
}

void parallel_write_plain_trimmed(const std::string &path,
                                  const RowDist &dist,
                                  const std::vector<std::uint8_t> &data,
                                  std::uint64_t original_size,
                                  std::uint64_t start_offset,
                                  bool atomic_io,
                                  MPI_Comm comm)
{
    const std::uint64_t base = (std::uint64_t)dist.offset_rows * dist.grid_size;
    const std::uint64_t len = (std::uint64_t)dist.local_rows * dist.grid_size;
    const std::uint64_t end = base + len;

    // Nur Bereich [start_offset, start_offset + original_size) in die Datei schreiben
    const std::uint64_t L = std::max(base, start_offset);
    const std::uint64_t R = std::min(end, start_offset + original_size);

    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(),
                            MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &fh),
              "open plain write", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity plain write", comm);

    if (R > L)
    {
        const std::uint64_t wr_bytes = R - L;
        const std::size_t src_off = (std::size_t)(L - base);
        const MPI_Offset file_off = (MPI_Offset)(L - start_offset);
        mpi_check(MPI_File_write_at(fh, file_off,
                                    data.data() + src_off,
                                    (int)wr_bytes, MPI_BYTE, MPI_STATUS_IGNORE),
                  "write_at plain(centered)", comm);
    }

    mpi_check(MPI_File_close(&fh), "close plain write", comm);

    // Dateigröße einmalig exakt setzen (kein Barrier nötig)
    int r0 = 0;
    MPI_Comm_rank(comm, &r0);
    if (r0 == 0)
    {
        MPI_File f2;
        mpi_check(MPI_File_open(MPI_COMM_SELF, path.c_str(),
                                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f2),
                  "open plain resize", MPI_COMM_SELF);
        mpi_check(MPI_File_set_size(f2, (MPI_Offset)original_size),
                  "set_size plain post", MPI_COMM_SELF);
        mpi_check(MPI_File_close(&f2), "close plain resize", MPI_COMM_SELF);
    }
}

void dump_frame_parallel(const std::string &path,
                         const RowDist &dist,
                         const std::vector<std::uint8_t> &data,
                         std::size_t padded_size_bytes,
                         bool atomic_io,
                         MPI_Comm comm)
{
    MPI_File fh;
    mpi_check(MPI_File_open(comm, path.c_str(),
                            MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &fh),
              "open frame", comm);
    mpi_check(MPI_File_set_atomicity(fh, atomic_io ? 1 : 0), "set_atomicity frame", comm);

    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
    {
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

bool read_meta_rank0(const std::string &meta_path,
                     std::uint64_t &original_size,
                     std::uint32_t &grid_size_n,
                     std::uint64_t &start_offset)
{
    std::ifstream in(meta_path, std::ios::binary);
    if (!in)
        return false;

    // Versuchen: magic + version. Wenn nicht vorhanden, alten zweifeldigen Header lesen.
    std::uint32_t magic = 0, version = 0;
    in.read(reinterpret_cast<char *>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char *>(&version), sizeof(version));
    if (!in)
        return false;

    if (magic == 0x48505031)
    { // "HPP1"
        // v1/v2
        in.read(reinterpret_cast<char *>(&original_size), sizeof(original_size));
        in.read(reinterpret_cast<char *>(&grid_size_n), sizeof(grid_size_n));
        if (!in)
            return false;
        if (version >= 2)
        {
            in.read(reinterpret_cast<char *>(&start_offset), sizeof(start_offset));
            if (!in)
                return false;
        }
        else
        {
            start_offset = 0;
        }
        return true;
    }
    else
    {
        // Rückfälliger Pfad: alte Meta ohne magic/version (nur 2 Felder)
        // Wir haben schon 8 Bytes (oder 4+4) verbraucht → zurückspulen:
        in.clear();
        in.seekg(0, std::ios::beg);
        in.read(reinterpret_cast<char *>(&original_size), sizeof(original_size));
        in.read(reinterpret_cast<char *>(&grid_size_n), sizeof(grid_size_n));
        if (!in)
            return false;
        start_offset = 0;
        return true;
    }
}

bool write_meta_rank0(const std::string &meta_path,
                      std::uint64_t original_size,
                      std::uint32_t grid_size_n,
                      std::uint64_t start_offset)
{
    std::ofstream out(meta_path, std::ios::binary | std::ios::trunc);
    if (!out)
        return false;

    // Write magic + versioned header (v2 includes start_offset)
    const std::uint32_t magic = 0x48505031; // "HPP1"
    const std::uint32_t version = 2;        // v2: with start_offset
    out.write(reinterpret_cast<const char *>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char *>(&version), sizeof(version));

    // Common fields
    out.write(reinterpret_cast<const char *>(&original_size), sizeof(original_size));
    out.write(reinterpret_cast<const char *>(&grid_size_n), sizeof(grid_size_n));
    // v2 field
    out.write(reinterpret_cast<const char *>(&start_offset), sizeof(start_offset));

    return static_cast<bool>(out);
}