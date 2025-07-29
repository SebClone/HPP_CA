# HPP Cellular Automaton Bit-State Documentation

This document provides a detailed explanation of the bit-level representation used in the HPP (Hardy–Pomeau–de Pazzis) cellular automaton.

## Bit-State Format

Each cell in the grid is represented by a `uint8_t` value:

```
uint8_t Bit = xxxxxxxx
```

- The first 4 bits (from the left) are the **high bits**.
- The last 4 bits are the **low bits**.

```
uint8_t Bit = xxxx|xxxx
                   ↑     ↑
                High   Low
```

### Wall Bit (High Bits)

- The last bit of the high nibble indicates whether a wall (obstacle) exists in the cell:
  - `0`: No wall
  - `1`: Wall present

**Examples**:

- No wall: `Bit = xxx0|xxxx`
- Wall: `Bit = xxx1|xxxx`

### Particle Directions (Low Bits)

The low bits represent the presence of particles in the four cardinal directions:

- N (North): 0b1000
- E (East): 0b0100
- S (South): 0b0010
- W (West): 0b0001

Each cell can hold at most one particle in each direction.

**Examples**:

- One particle:

  - `Bit = xxxx|1000` → North
  - `Bit = xxxx|0100` → East
  - `Bit = xxxx|0010` → South
  - `Bit = xxxx|0001` → West

- Two particles:

  - `Bit = xxxx|1001` → North-West
  - `Bit = xxxx|1010` → North-South
  - `Bit = xxxx|1100` → North-East
  - `Bit = xxxx|0101` → East-West
  - `Bit = xxxx|0110` → East-South
  - `Bit = xxxx|0011` → South-West

- Three particles:

  - `Bit = xxxx|1011` → N-S-W
  - `Bit = xxxx|1101` → N-E-W
  - `Bit = xxxx|1110` → N-E-S
  - `Bit = xxxx|0111` → E-S-W

- Four particles:
  - `Bit = xxxx|1111` → N-E-S-W

> Each element in the NxN automaton grid corresponds to a cell encoded in this format.

## Collision Rules

- When two particles collide head-on (e.g., W → ← E), they transform into a vertical pair (N + S), and vice versa.

**Examples**:

- `Bit = xxx0|1010` (N + S) ⟶ `Bit = xxx0|0101` (E + W)
- `Bit = xxx0|0101` (E + W) ⟶ `Bit = xxx0|1010` (N + S)

- Particles do **not** collide in other configurations. For 3 or 4 particles, the state remains unchanged.

> **Note**: Cells with the wall bit set (`xxx1|xxxx`) do **not** participate in collisions and follow a separate transition logic.

## Propagation Rules

After the collision step (or before—it is symmetric due to cyclic behavior), particles propagate to their neighboring cells in the direction they are moving.

- N particle → upper neighbor
- E particle → right neighbor
- S particle → lower neighbor
- W particle → left neighbor

**Example Implementation**:

```cpp
if (current_cell & 0b00001000) // North particle present
    upper_neighbor |= 0b00001000; // Add North particle to upper neighbor
```

> This logic ensures particle states are added to neighbors regardless of existing directions.

## Reflection Rules

Cells with obstacles reflect particles by reversing their direction (180° rotation):

- N → S
- E → W
- S → N
- W → E

**Examples**:

- `Bit = xxx1|1000` (N) ⟶ `Bit = xxx1|0010` (S)
- `Bit = xxx1|0100` (E) ⟶ `Bit = xxx1|0001` (W)
- `Bit = xxx1|0010` (S) ⟶ `Bit = xxx1|1000` (N)
- `Bit = xxx1|0001` (W) ⟶ `Bit = xxx1|0100` (E)

> Reflections are handled separately from collisions due to differing physical behavior.
