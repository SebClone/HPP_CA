# %% Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib import colors
import seaborn as sns
from pathlib import Path

saveAnimation = True  # Set to False to disable saving animations

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

def animate_frames(frames_dict, gridsize, outfile='grid_frames_visualization.mp4', fps=20, cmap='gray'):
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
def animate_frame_diffs(frames_dict, gridsize, outfile='frames_diff_visualization.mp4', fps=20, cmap='bwr'):
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
grid_size = int(np.sqrt(example_frame.size))
print(f"Grid size: {grid_size}x{grid_size}")

# %% Print first and last frame as text
print("\n--- First Frame Text ---\n")
first_frame = frames['frame_000000']
print(binary_to_text(first_frame))

print("\n--- Last Frame Text ---\n")
last_frame = frames['frame_000999']
print(binary_to_text(last_frame))


# %% Animation of frames
frame_0 = np.fromfile("frames/frame_000000.bin", dtype=np.uint8)
plot_frame(frame_0, "Frame 0", grid_size)
frame_999 = np.fromfile("frames/frame_000999.bin", dtype=np.uint8)
plot_frame(frame_999, "Frame 999", grid_size)

# %%
# Animate all frames and save as MP4
if saveAnimation:
    print("Saving animations...")
    animate_frames(frames, grid_size, outfile='visualization_results/grid_frames_visualization.mp4', fps=20, cmap='gray')
    animate_frame_diffs(frames, grid_size, outfile='visualization_results/frames_diff_visualization.mp4', fps=20, cmap='bwr')

# %% Plot difference between first and last frame
frame0 = frames['frame_000000']
frame999 = frames['frame_000999']
plot_frame_diff(frame0, frame999, grid_size, title="Pixelweise Differenz (uint8)")
