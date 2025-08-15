# %% Imports
import numpy as np
import pandas as pd
import io
from PIL import Image, UnidentifiedImageError
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib import colors
import seaborn as sns
from pathlib import Path
from PIL import Image
import struct
from typing import Optional

# --- Read encrypted_full_meta and parse original size and grid size ---
meta_path = Path("encrypted_full.meta")
with open(meta_path, "rb") as f:
    original_size_bytes = f.read(8)
    original_size = struct.unpack("<Q", original_size_bytes)[0]
    Nmeta_bytes = f.read(4)
    Nmeta = struct.unpack("<I", Nmeta_bytes)[0]
print(f"Original size from encrypted_full_meta: {original_size}")
grid_size = int(Nmeta)
print(f"Grid size from encrypted_full_meta: {grid_size}")

isText = False  # Set to True to print frames as text
saveAnimation = False  # Set to False to disable saving animations

# %% Functions
def plot_benchmark_data(processes, times):
    sns.set(style="whitegrid", context="talk")
    plt.figure(figsize=(10, 6))
    sns.lineplot(x=processes, y=times, marker="o", color="b")
    plt.title('Benchmark Results', fontsize=22)
    plt.xlabel('Number of Processes', fontsize=18)
    plt.ylabel('Time [s]', fontsize=18)
    plt.xticks(processes, fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)
    plt.savefig('visualization_results/benchmark_results.png', dpi=300, bbox_inches='tight')
    plt.show()

def binary_to_text(binary_data):
    # Convert binary data to bytes
    raw_data = binary_data.tobytes()
    # Remove common padding (NUL bytes) and try UTF-8, then fall back
    raw_data = raw_data.replace(b"\x00", b"")
    try:
        return raw_data.decode("utf-8")
    except UnicodeDecodeError:
        return raw_data.decode("latin-1", errors="replace")
    

def plot_frame(frame_data, title, gridsize, cmap='gray'):
    arr_2d = frame_data.reshape((gridsize, gridsize))
    plt.imshow(arr_2d.astype(np.int16), cmap=cmap, interpolation='nearest', vmin=0, vmax=255)
    plt.title(title)
    plt.colorbar(label="Value")
    plt.savefig(f"visualization_results/{title.replace(' ', '_')}.png", dpi=300, bbox_inches='tight')
    plt.show()

def grid_to_bytes_from_nibbles(cells_uint8: np.ndarray,
                               order: str = "hi-lo",
                               original_size: Optional[int] = None) -> bytes:
    """
    Rekonstruiert die ursprünglichen Datei-Bytes aus CA-Zellwerten (0..15 pro Zelle).
    Erwartet: 1 Byte pro Zelle im Dump (Wertbereich 0..15).
    order: 'hi-lo'  -> Byte = (cells[0]<<4) | cells[1]
           'lo-hi'  -> Byte = (cells[1]<<4) | cells[0]
    """
    # Nur die unteren 4 Bit verwenden – Sicherheitshalber maskieren
    cells = (cells_uint8 & 0x0F).astype(np.uint8)

    # Paare bilden: (0,1), (2,3), ...
    hi = cells[0::2]
    lo = cells[1::2]

    if order == "hi-lo":
        out = (hi << 4) | lo
    elif order == "lo-hi":
        out = (lo << 4) | hi
    else:
        raise ValueError("order must be 'hi-lo' or 'lo-hi'")

    if original_size is not None:
        out = out[:original_size]

    return out.tobytes()

def show_frame_decoded_from_grid(frame_cells: np.ndarray,
                                 N: int,
                                 original_size: int,
                                 order: str = "hi-lo",
                                 title: str = "decoded from grid"):
    """
    Baut aus einem Gitter-Frame (Zellenwerte) wieder die Datei-Bytes,
    decodiert sie als Bild und zeigt das RGB-Bild an.
    Hinweis: Das ergibt *nur* beim letzten Frame eines erfolgreichen Decrypts
    einen gültigen PNG/JPG-Stream.
    """
    # First try to treat frame_cells[:original_size] as complete image file
    try:
        file_bytes = frame_cells[:original_size].tobytes()
        im = Image.open(io.BytesIO(file_bytes)).convert("RGB")
        arr = np.array(im)
        plt.imshow(arr)
        plt.title(f"{title} ({im.width}×{im.height})")
        plt.axis("off")
        plt.show()
        return
    except UnidentifiedImageError:
        pass

    # Falls der Dump gepackt war (N*N/2 Bytes): erst Nibbles entpacken
    if frame_cells.size == (N * N) // 2:
        # hier liegen bereits gepaarte Nibbles je Byte vor -> direkt kürzen
        file_bytes = frame_cells[:original_size].tobytes()
    elif frame_cells.size == (N * N):
        # 1 Zelle pro Byte (Wert 0..15): zu Bytes packen
        file_bytes = grid_to_bytes_from_nibbles(frame_cells, order, original_size)
    else:
        raise ValueError(f"Unerwartete Framegröße: {frame_cells.size} (N={N})")

    # Bild aus Bytes decodieren
    im = Image.open(io.BytesIO(file_bytes)).convert("RGB")
    arr = np.array(im)
    plt.imshow(arr)
    plt.title(f"{title} ({im.width}×{im.height})")
    plt.axis("off")
    plt.show()

def plot_frame_diff(frame_a, frame_b, gridsize, title=None, cmap='bwr'):
    arr_a = frame_a.reshape((gridsize, gridsize)).astype(np.int16)
    arr_b = frame_b.reshape((gridsize, gridsize)).astype(np.int16)
    diff = arr_a - arr_b
    max_abs = int(np.max(np.abs(diff)))
    plt.imshow(
        diff,
        cmap=cmap,
        interpolation='nearest',
        vmin=-max_abs if max_abs > 0 else -1,
        vmax=max_abs if max_abs > 0 else 1
    )
    if title is None:
        plt.title("Differenz der Frames")
        filename = 'visualization_results/frame_diff.png'
    else:
        plt.title(title)
        filename = f"visualization_results/{title.replace(' ', '_')}.png"
    plt.colorbar(label="Delta")
    print("equal:", np.array_equal(frame_a.reshape((24,24)), frame_b.reshape((24,24))))
    print("max abs diff:", max_abs)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()


# Animation function for frames

def animate_frames(frames_dict, gridsize, outfile='visualization_results/grid_frames_visualization.mp4', fps=20, cmap='gray'):
    """
    Animates a sequence of frames stored in frames_dict and saves as MP4.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, FFMpegWriter

    keys = list(frames_dict.keys())
    first_frame = frames_dict[keys[0]].reshape((gridsize, gridsize)).astype(np.int16)
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(first_frame, cmap=cmap, interpolation='nearest', vmin=0, vmax=255)
    cbar = fig.colorbar(im, ax=ax, label="Value")
    title = ax.set_title(keys[0])

    def update(i):
        frame_name = keys[i]
        frame_arr = frames_dict[frame_name].reshape((gridsize, gridsize)).astype(np.int16)
        im.set_data(frame_arr)
        title.set_text(frame_name)
        return [im, title]

    anim = FuncAnimation(fig, update, frames=len(keys), interval=1000/fps, blit=False)
    writer = FFMpegWriter(fps=fps, codec='libx264')
    anim.save(outfile, writer=writer)
    plt.close(fig)
    print(f"Animation saved as {outfile}")


# Animation of frame differences
def animate_frame_diffs(frames_dict, gridsize, outfile='visualization_results/frames_diff_visualization.mp4', fps=20, cmap='bwr'):
    """
    Animates the differences between consecutive frames and saves as MP4.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, FFMpegWriter

    keys = list(frames_dict.keys())
    if len(keys) < 2:
        raise ValueError("Need at least two frames to compute differences")
    # Compute diffs
    diffs = []
    for i in range(1, len(keys)):
        arr_prev = frames_dict[keys[i-1]].reshape((gridsize, gridsize)).astype(np.int16)
        arr_i = frames_dict[keys[i]].reshape((gridsize, gridsize)).astype(np.int16)
        diffs.append(arr_i - arr_prev)
    max_abs = int(np.max([np.max(np.abs(d)) for d in diffs]))
    vmin = -max_abs if max_abs > 0 else -1
    vmax = max_abs if max_abs > 0 else 1

    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(diffs[0], cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax, label="Delta")
    title = ax.set_title(f"Δ {keys[1]} − {keys[0]}")

    def update(i):
        im.set_data(diffs[i])
        title.set_text(f"Δ {keys[i+1]} − {keys[i]}")
        return [im, title]

    anim = FuncAnimation(fig, update, frames=len(diffs), interval=1000/fps, blit=False)
    writer = FFMpegWriter(fps=fps, codec='libx264')
    anim.save(outfile, writer=writer)
    plt.close(fig)
    print(f"Animation saved as {outfile}")

# %% Benchmark data
# Load benchmarking results
benchmark_results = pd.read_csv(
    'hpp_benchmarking_results.txt',
    sep='|',
    skipinitialspace=True
)
benchmark_results.columns = benchmark_results.columns.str.strip()
processes = benchmark_results['Prozesse']
times = benchmark_results['Laufzeit (s)']
plot_benchmark_data(processes, times)


# %% Load grid data
frames = {}
from pathlib import Path
folder_path = Path("frames")
# Find all .bin files, sort them, and load into frames dict
bin_files = sorted(folder_path.glob("*.bin"))
for bin_file in bin_files:
    # Remove suffix to get frame name
    frame_name = bin_file.stem
    arr = np.fromfile(bin_file, dtype=np.uint8)
    frames[frame_name] = arr

# %% Calculate grid size
example_frame = next(iter(frames.values()))
# Read image_message.png to get width and height
with Image.open("image_message.png") as img:
    width, height = img.size
# grid_size is now from meta, so no need to override here
picture_size = [width, height]
print(f"Grid size: {grid_size}x{grid_size}")
print(f"Picture size: {picture_size[0]}x{picture_size[1]}")


# %% Print first and last frame
first_frame = frames['frame_000000']
last_frame = frames['frame_000999']
if isText:
    print("\n--- First Frame Text ---\n")
    print(binary_to_text(first_frame))

    print("\n--- Last Frame Text ---\n")
    print(binary_to_text(last_frame))

    frame_0 = np.fromfile("frames/frame_000000.bin", dtype=np.uint8)
    plot_frame(frame_0, "Frame 0", grid_size)
    frame_999 = np.fromfile("frames/frame_000999.bin", dtype=np.uint8)
    plot_frame(frame_999, "Frame 999", grid_size)
    frame0 = frames['frame_000000']
    frame999 = frames['frame_000999']
    plot_frame_diff(frame0, frame999, grid_size, title="Pixelweise Differenz (uint8)")
else:
    print("Data is an image")
    show_frame_decoded_from_grid(last_frame, N=grid_size,
                                original_size=original_size,
                                order="hi-lo",
                                title="Frame 999 (decoded)")

# Animate all frames and save as MP4
if saveAnimation & isText:
    print("Saving animations...")
    animate_frames(frames, grid_size, outfile='visualization_results/grid_frames_visualization.mp4', fps=20, cmap='gray')
    animate_frame_diffs(frames, grid_size, outfile='visualization_results/frames_diff_visualization.mp4', fps=20, cmap='bwr')


