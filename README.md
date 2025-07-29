# HPP_CA

This project is part of the "Parallel Computing" lecture in the summer semester of 2025.

## Project Overview

This implementation uses a HPP (Hardy–Pomeau–de Pazzis) cellular automaton to encrypt and decrypt text messages with toroidal boundary conditions.

### File Descriptions

#### `HPP_lib.hpp`
Declares all functions necessary for encryption and decryption using the HPP cellular automaton.

#### `HPP_lib.cpp`
Defines the functions declared in `HPP_lib.hpp`, providing the core logic used in the encryption and decryption algorithms.

#### `encryption_HPP.cpp`
Implements the algorithm for encrypting a text message using the HPP cellular automaton. The steps include:
- Generating a wall bitmask using a random method, assigning wall bits to approximately 10% of the `uint8_t` elements.
- Exporting the wall bitmask as a `.key` file for later use during decryption.
- Loading the original message from a `.txt` file, converting it into a binary format, and reshaping it into a 2D matrix.
- Running the HPP algorithm for **N** steps:
  - Particle collision
  - Particle propagation
  - Particle reflection on walls
- Saving the encrypted message as a `.txt` file, ready for transmission to the recipient.

#### `decryption_HPP.cpp`
Implements the algorithm for decrypting a text message using the HPP cellular automaton. The steps include:
- Loading the wall bitmask from the previously exported `.key` file.
- Loading the encrypted message from a `.txt` file, converting it into binary format, and reshaping it into a 2D matrix.
- Running the HPP algorithm in reverse for **N** steps:
  - Particle reflection on walls
  - Particle propagation (in reverse)
  - Particle collision
- Saving the decrypted message as a `.txt` file.

## Notes
- Toroidal boundary conditions ensure the grid wraps around at the edges.
- The encryption is reversible only if the same wall bitmask and number of iterations **N** are used.
