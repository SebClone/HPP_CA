__(englisch Version below)__

# Dokumentation für das Projekt "HPP-Automat" im Modul "Parallel Computing"

Von Samuel Orth und Sebastian Roth.

## Tldr:
Symmetrischer Verschlüsselungsalgorithmus auf Basis eines zellulären Gasautomatens.
Die Parallelisierung erfolgt mittels MPI und domain-decomposition in 2D Blöcken mit Halo-Zellen Kommunikation der Nachbarn.
Nachbarschaftskommunikation und Berechnung des unabhängigen inneren Kerns erfolgen überlappend.
Für den inneren Kern wurde zusätzlich OpenMP verwendet.

## Inhaltverzeichnis:
0. Projektinfo
1. Nutzung
2. Projektstuktur
3. Umsetzung des HPP-Automaten
4. Designentscheidungen
5. Parallelisierung
6. Allgemeine Hinweise zum Code

## 0) Projektinfo
Dieses repository implementiert einen Lattice Gas Cellular Automaton - ein symmetrischer Verschlüssenlungsalgorithmus nach dem Vorbild von 
Laurent Signac's ["Lattice Gas Symmetric Cryptography"](https://doi.org/10.48550/arXiv.1306.1519).
Nachrichten im Binärformatwerden in einer Matrix aus uint8_t Variablen abgespeichert und unter Einsatz spezieller Regeln (siehe 3.) iterativ verändert - die Bits der Nachricht bewegen sich wie "Teilchen" in einer Gitterbox.
Als Schlüssel fungiert die Position sogenannter Wandzellen im Gitter, die nur der Empfänger der Nachricht erhält und über "rückwärts abspulen" des Verschlüsselten Zustandes die Ausgangsnachricht entschlüsseln kann.
Ohne genaue Position der Wandzellen ist keine eindeutige Rekonstruktion der Nachricht möglich.

## 1) Nutzung
Wir nutzen ein Makefile (/Makefile) und eine Konfigurationsdatei (/include/config.hpp, /src/config.cpp).
Einstellbar in der Config sind:
- iterations             (typ: integer , default: 1000,                         Aufgabe: Anzahl der Iterationen des Haupt-loops)
- grid_size              (typ: integer , default: automatisch,                  Aufgabe: Größe des Gitters, in der die Nachricht eingelesen wird)
- wall_density           (typ: double  , default: 0.10,                         Aufgabe: Verhältnis Wandzellen zu freien Gitterzellen (Bereich [0.0, 1.0]))
- seed                   (typ: uint64_t, default: random,                       Aufgabe: Seed für die randomisierte Verteilung der Wanzellen)
- dump_frames            (typ: boolean , default: false,                        Aufgabe: Speichert zwischendurch "Snapshots" des aktuellen Grids)
- frame_interval         (typ: integer , default: 10,                           Aufgabe: Abstand in dem Snapshots gespeichert werden)
- input                  (typ: string  , default: "data/message.txt",           Aufgabe: Pfad zur Originalnachricht)
- enc_bin                (typ: string  , default: "data/encrypted_full.bin",    Aufgabe: Pfad zur verschlüsselten Nachricht)
- meta                   (typ: string  , default: "data/encrypted_full.meta",   Aufgabe: Pfad zu den Metadaten)
- key                    (typ: string  , default: "data/wall_mask.key",         Aufgabe: Pfad zum Key)
- output                 (typ: string  , default: "data/decrypted_message.txt", Aufgabe: Pfad zur Entschlüsselten Nachricht)

Einstellungen im Makefile:
- mpic++ Compiler
- build target heißt "hpp_mpi_app"
- Basis-Flags: -std=c++20 -Wall -Wextra -Wpedantic
- Optimiert standardmäßig it -O3 (es sei denn BUILD=debug ist gesetzt)
- OpenMP wird genutzt (es sei denn OMP=0 wird gesetzt)
- Laufzeitparameter: - NP=4 (Zahl MPI-Prozesse)
                     - RUN_ARGS=("") geht an mpirun
                     - APP_ARGS=("") geht ans Programm
- CLI Flag: --encrypt und --decrypt legen Modus fest (Übergabe an APP_ARGS oder einfach make encrypt bzw. make decrypt)

Beispielverwendung:
- make encrypt
- make run NP=8 RUN_ARGS='-x OMP_NUM_THREADS=4' APP_ARGS='--decrypt"

## 2) Projektstuktur
* /include/ enthält alle .hpp Header.
* /src/ enthält alle .cpp Quellfiles.
* /.build/ enthält alle automatisch erzeugten build-objekte.
* /results/ enthält alle von uns produzierten ERgebnisse, wie Graphen, Tabellen etc.
* /data/ enthält alle input und output files, die wichtig für den Encryption und Decryption Modus sind. 

## 3) Umsetzung des HPP-Automaten
Die Implementation aller HPP-Regeln befindet sich in /hpp_rules.hpp und /hpp_rules.cpp.
Der HPP-Cellular Automata ist ein Zellautomat mit einer van Neumann NAchbarschaft. Wir haben uns für eine Torus-charakteristik entschieden. Somit gibt es keine Enden des Gitters und die Randspalten/-zeilen sind in direkter Nachbarschaft zueinander.
Der Automat besteht aus einer 2D-Matrix (dtype=uint8_t) der größe grid_size x grid_size, die per Termminal eingabe definiert wird.

Eine belibige Nachricht/ Datei wird als bytes/ bits eingelesen und in die 2D_Matrix, also das grid, geschrieben. 
Die Nachricht kann mit einem start_offset in das grid geschrieben werden.

Für den HPP-Zellautomaten sind nur die low-bits relevant. Diese sind wie folgt codiert:
- One particle:

  - `Bit = xxxx|1000` → Nord
  - `Bit = xxxx|0100` → Ost
  - `Bit = xxxx|0010` → Süd
  - `Bit = xxxx|0001` → West

Für die Simmulation werden die HPP-Regeln verwendet. Diese sind in [HPP States](HPP_states.md) detailiert beschrieben.
Die Wall-Bit-Maske wird dabei parallel als eine 2D-matrix (dtype=bolean) der selben größe angelegt. Durch die seperate Behandlung der Mask wird verhindert, dass Informationen der original Datei verloren gehen , wenn das erste high-bit manipuliert wird.

Zur Encryption wird Collision -> Propagation -> Refelction nacheinander auf das grid angewendet.
Zur Decryption wird Reflection -> Inverse-Propagation -> Collision nacheinander auf das grid angewendet.

***HPP-Regeln kurz***
- Collision: Zwei frontal aufeinander treffende particle (Nord-Süd/ Ost-West) werden um 90 grad rotiert 
    - `Bit = xxx0|1010` (N + S) ⟶ `Bit = xxx0|0101` (E + W)
- Reflection: Ist in der Zelle eine Wand wird das Partikel um 180 grad gedreht (reflektiert)
    - `Bit = xxx1|1000` (N) ⟶ `Bit = xxx1|0010` (S)
- Propagation: Ein bit wird entlang seiner Richtung an die entsprechende NAchbarzelle weitergegeben
    - Nord      -> Oben
    - Ost       -> Rechts
    - Süd       -> Unten
    - West      -> Links
- Inverse-Propagation: Ein bit aus der jeweiligen Nachbarzelle wird entgegen seiner Richtung "zurückgegeben".
    - Oben      -> Nord
    - Rechts    -> Ost
    - Unten     -> Süd
    - Links     -> West

## 4) Parallelisierung
Wir verwenden MPI zur Zerlegung des Grids in 2D-Blöcke (domain-decomposition).
Aufgrund der Torus-Topologie des Grids besitzt jeder Block 4 Nachbarn (up/down und right/left).
Für solche Aufgaben stellt MPI die MPI_Cart Befehlsgruppe zur Verfügung.
Für den Übergang von "Teilchen" zwischen Domänen müssen diese miteinander kommunizieren.
Dies geschieht über das Einführen von zusätzlichen Ghost-Zellen (auch Halo-Zellen) an den Rändern jedes 2D-Blocks, in die in jedem Iterationsschritt die Werte den vier Nachbarn neu eingespielt werden (MPI_Isend, MPI_Irecv jeder Randspalte/Randzeile).
Mit der Nutzung von MPI_Isend und MPI_Irecv kann die Kommunikation der Halo-Zellen überlappend mit der Berechnung der inneren (von den Halo-Zellen unabhängigen) Zellen erfolgen.
Die Synchronisation erfolgt mittels eines nachgeschalteten MPI_Waitall.
Nach Erhalt der Informationen aus benachbarten Randzellen werden können die Randzellen jedes Blockes per "pragma omp simd" berechnet werden.
Alle Anwendung der HPP-Regeln erfolgt im Algorithmus für jeden Rank Zellenweise, jedoch müssen sich die Bits so Verhalten, als wenn sie gleichzeitig propagiert würden.
Daher wird mit einem Doppel-Buffer System gearbeitet (Berechnung auf Basis des aktuellen Zustandes in Buffer A, schreiben des neuen Zustandes in Buffer B --> wechsel der Buffer nach jeder Iteration).

Fußnote: MPI wurde auch zur Parallelsierung für einlesen und schreiben der Daten in 1D-Zeilenstreifen verwendet (siehe io.cpp)
Jedoch ist bei standard Gridgrößen und Nachrichtenlängen die Laufzeit zu >98% im Haupt-Loop, das Hauptaugemnerk liegt hier also auf der Parallelisierung dieses Loops.

## 5) Designentscheidungen
- Wir speichern zusätzliche Metadaten (Originalgröße der Nachricht, Beginn der Nachricht im Grid, Gridgröße) in einem separaten File. Denkbar wäre auch die Verwendung eines "Stopp-bytes" um Ende (und ggf. Beginn) der Nachricht zu markieren.
- Die Gridgröße richtet sich standarmäßig nach der Größe der Nachricht. Feste Gridgrößen sind jedoch ebenfalls einstellbar und wurden zur Erzeugung der Ergebnisse verwendet.
- Wir verwenden Standardmäßig Textnachrichten zum testen der Verschlüsselung. Da der Algorithmus mit Binärdaten arbeitet ließen sich im Prinzipa auch andere Datenformate verschlüsseln.

## 6) Allgemeine Hinweise zum Code
-Die include Befehle jedes Files erfolgen nach dem "include what you use" Prinzip.
-Das etwas längliche umrechnen und mappen von 1D <-> 2D mit MPI_Alltoallw ist dem geschuldet, dass die Funktionen in io.cpp zuest für 1D Zeilenstreifen ausgelegt waren. Somit war ein Adapter nötig.




## English HPP
The implementation of all HPP rules is located in /hpp_rules.hpp and /hpp_rules.cpp.
The HPP cellular automata is a cellular automata with a von Neumann neighborhood. We opted for a torus characteristic. This means that there are no ends to the grid and the edge columns/rows are in direct proximity to each other.
The automaton is a 2D matrix (dtype=uint8_t) of size grid_size x grid_size, which is defined by terminal input.

Any message/file is read in as bytes/bits and written to the 2D matrix, i.e., the grid. 
The message can be written to the grid with a start_offset.

Only the low bits are relevant for the HPP cellular automaton. These are encoded as follows:
  - `Bit = xxxx|1000` → North
  - `Bit = xxxx|0100` → East
  - `Bit = xxxx|0010` → South
  - `Bit = xxxx|0001` → West

The HPP rules are used for the simulation. These are described in detail in [HPP States](HPP_states.md).
The wall bit mask is created in parallel as a 2D matrix (dtype=boolean) of the same size. Separate handling of the mask prevents information from the original file from being lost when the first high bit is manipulated.

For encryption, collision -> propagation -> reflection is applied to the grid in succession.
For decryption, reflection -> inverse propagation -> collision is applied to the grid in succession.

***HPP-ruels brief***
- Collision: Two particles colliding head-on (north-south/east-west) are rotated by 90 degrees
    - `Bit = xxx0|1010` (N + S) ⟶ `Bit = xxx0|0101` (E + W)
- Reflection: If there is a wall in the cell, the particle is rotated 180 degrees (reflected)
    - `Bit = xxx1|1000` (N) ⟶ `Bit = xxx1|0010` (S)
- Propagation:  A bit is passed along its direction to the corresponding neighboring cell
    - North      -> Up
    - East       -> Right
    - South       -> Down
    - West      -> Left
- Inverse-Propagation: A bit from the respective neighboring cell is “returned” in the opposite direction.
    - Up      -> North
    - Right    -> East
    - Down     -> South
    - Left     -> West




