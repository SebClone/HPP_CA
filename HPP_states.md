Bit-State explanation:
uint8_t Bit = xxxxxxxx - The first 4 bits are the so called high bits and the last 4 bits are called low bits
e.g. uint8_t Bit = xxxx|xxxx
High|low - The the last bit of the High bits determin wehter there is an wall/obstacle in that cell or not
e.g. NO-wall uint8_t Bit = xxx0|xxxx
High|low
e.g. wall uint8_t Bit = xxx1|xxxx
High|low - The low bits determin wehter a partical is in an certain direction
There are 4 directions North (N), East (E), South (S) and West (W). Per cell only one partical in each direction is allowed. So a Cell can have up to 1 partical in (N),(E),(S) and (W). The bits are then defined as follows.
e.g. uint8_t Bit = xxxx|NESW
High|low
One Partical:
e.g. uint8_t Bit = xxxx|1000 North  
 e.g. uint8_t Bit = xxxx|0100 East
e.g. uint8_t Bit = xxxx|0010 South
e.g. uint8_t Bit = xxxx|0001 West
Two Paticals:
e.g. uint8_t Bit = xxxx|1001 North-West
e.g. uint8_t Bit = xxxx|1010 North-South
e.g. uint8_t Bit = xxxx|1100 Noth-East
e.g. uint8_t Bit = xxxx|0101 East-West
e.g. uint8_t Bit = xxxx|0110 East-South
e.g. uint8_t Bit = xxxx|0011 South-West
Three Particals:
e.g. uint8_t Bit = xxxx|1011 North-South-West
e.g. uint8_t Bit = xxxx|1101 North-East-West
e.g. uint8_t Bit = xxxx|1110 North-East-South
e.g. uint8_t Bit = xxxx|0111 East-South-West
Four Particals:
e.g. uint8_t Bit = xxxx|1111 North-East-South-West > The cellular automata is an NxN-matrix where each element of the matrix corresponds to an cell that looks like this e.g. uint8_t Bit = xxxx|1111 North-East-South-West. So each cell can have an obstacle/wall and one of each partical type

Collision: - When two particals collide frontal so (W)-> <-(E) then they change their directions to (N) and (S)
The state transition would look like this:
e.g. uint8_t Bit = xxx0|1010 North-South ------> e.g. uint8_t Bit = xxx0|0101 East-West
e.g. uint8_t Bit = xxx0|0101 East-West ------> e.g. uint8_t Bit = xxx0|1010 North-South - Other then this partical will not collide. So also if there are three particals the state will remain the same.
When four particals are in one cell theoretically they will collide but agian into a North-East-South-West configuration
! Cells with an wall bit set to 1, do not take part in the collison step. Since they follow an diffrent state transition table. Also physically the particals never collide with each other only with the wall

Propagation: - After or bevore, it does not matter, the particals collide they will propagate to the next neighbouring cell.
The difrence in the order comes down to the assumtion of the starting point since it is an cyclic iteration.
Each partical moves to the direction they are in, so an South partical will move to the bottom neighbour and so forth. > It will be implementet like this:
if (current_cell & 0b00001000) (Checks wehter in the current cell is an north partical)
upper neighbour =| 0b00001000 (Gives the upper neighbour the north partical. This functions as an addition. So the upper neighbour gets an north particel independent wether it has an east, south or west partical)

Reflection: - In this lattice there are obstical or walls, as touchted above. on these a partical is not only deflected but relfected. So it chages its direction 180˚.
Our method lets partical propagate an cell that contains an obstacle/wall. Yet it is not an colission because its state transition is diffrent. Hence it makes sence to implement a diffrent function that only handels reflections.
The state transitions would look like this:
e.g. uint8_t Bit = xxx1|1000 North ------> e.g. uint8_t Bit = xxx1|0010 South
e.g. uint8_t Bit = xxx1|0100 East ------> e.g. uint8_t Bit = xxx1|0001 West
e.g. uint8_t Bit = xxx1|0010 South ------> e.g. uint8_t Bit = xxx1|1000 North
e.g. uint8_t Bit = xxx1|0001 West ------> e.g. uint8_t Bit = xxx1|0100 East

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
