from pathlib import Path
import struct
import numpy as np
from PIL import Image
import imageio.v2 as imageio  # pip install imageio

# --- Meta lesen (original_size, N) ---
with open("encrypted_full.meta", "rb") as f:
    original_size = struct.unpack("<Q", f.read(8))[0]
    N = struct.unpack("<I", f.read(4))[0]
print(f"original_size={original_size}, N={N}")

in_dir  = Path("frames/binary")
out_dir = Path("frames/png"); out_dir.mkdir(parents=True, exist_ok=True)

png_paths = []
for bin_path in sorted(in_dir.glob("*.bin")):
    buf = np.fromfile(bin_path, dtype=np.uint8)

    # Erwartete Länge für N×N:
    needed = N * N
    if buf.size < needed:
        # auffüllen (wie PGM‑Variante) – sonst hätten wir keine volle Matrix
        buf = np.pad(buf, (0, needed - buf.size), constant_values=0)
    elif buf.size > needed:
        # falls mehr Bytes drin sind (z.B. N gerundet): abschneiden
        buf = buf[:needed]

    arr = buf.reshape(N, N)

    # Als **gültige PNG** in Graustufe („L“) speichern
    im = Image.fromarray(arr, mode="L")
    out_path = out_dir / (bin_path.stem + ".png")
    im.save(out_path, format="PNG")
    png_paths.append(out_path)
    print("wrote", out_path)


# --- optional: Animation als MP4 (streamend, speicherschonend) ---
make_mp4 = False
if make_mp4 and png_paths:
    mp4_path = out_dir / "frames.mp4"
    try:
        import imageio.v2 as imageio
        with imageio.get_writer(
            mp4_path,
            format="FFMPEG",
            mode="I",
            fps=20,
            codec="libx264",
            output_params=["-pix_fmt", "yuv420p"],  # breit kompatibel
        ) as writer:
            for p in png_paths:
                frame = imageio.imread(p)  # nur EIN Frame im RAM
                writer.append_data(frame)
        print("animation:", mp4_path)
    except Exception as e:
        print("Konnte MP4 nicht schreiben:", e)
        print("Tipp: `pip install imageio-ffmpeg` oder alternativ als GIF speichern.")