__(englisch Version below)__

# Dokumentation für das Projekt "HPP-Automaton" im Modul parallel computing.
Von Samuel Orth und Sebastian Roth.

## TldR: 

## Inhaltverzeichnis:
0. Projektinfo
1. Nutzung
2. Projektstuktur
3. Umsetzung des HPP-Automaten
4. Desing-Choices
5. Parallelisierung
6. Allgemeine Hinweise zum Code

## 0)
Dieses repository implementiert einen Lattice Gas Cellular Automaton - ein symmetrischer (...) nach dem Vorbild von (...).
Im Rahmen dieses Projekt wurde die encryption/decryption mithilfe von MPI und OpenMP parallelisisert

## 1) Nutzung
Wir nutzen ein Makefile (/Makefile) und eine Konfigurationsdatei (/include/config.hpp, /src/config.cpp).
Einstellbar in der Config sind:
- iterations             (typ: integer , default: 1000, Aufgabe: Anzahl der Iterationen des Haupt-loops)
- grid_size              (typ: integer , default: automatisch, Aufgabe: Größe des Gitters, in der die Nachricht eingelesen wird)
- wall_density           (typ: double  , default: 0.10, Aufgabe: Verhältnis Wandzellen zu freien Gitterzellen (Bereich [0.0, 1.0]))
- seed                   (typ: uint64_t, default: random, Aufgabe: Seed für die randomisierte Verteilung der Wanzellen)
- dump_frames            (typ: boolean , default: false, Aufgabe: Speichert zwischendurch "Snapshots" des aktuellen Grids)
- frame_interval         (typ: integer , default: 10, Aufgabe: Abstand in dem Snapshots gespeichert werden)
- input                  (typ: string  , default: "", Aufgabe: Pfad zur Originalnachricht (für Encryption))
- enc_bin                (typ: string  , default:, Aufgabe: Pfad zur verschlüsselten Nachricht (für Encryption und Decryption))
- meta                   (typ: string  , default:, Aufgabe: Pfad zu den Metadaten (für Encryption und Decryption))
- key                    (typ: string  , default:, Aufgabe: Pfad zum key (für Encryption und Decryption))
- output                 (typ: string  , default:, Aufgabe: Pfad zur Entschlüsselten Nachricht (für Decryption))

Einstellungen im Makefile:
-mpic++ Compiler
-build target heißt "hpp_mpi_app"
-basis flags für (...) -std=c++20 -Wall -Wextra -Wpedantic
-OpenMP:
-Laufzeitparameter: 
-CLI Flag:

Beispielverwendung:
make encrypt
make run NP=8 RUN_ARGS='-x OMP_NUM_THREADS=4' APP_ARGS='--decrypt"



## 2) Projektstuktur
* /include/ enthält alle .hpp Header.
* /src/ enthält alle .cpp Quellfiles
* /.build/ enthält alle automatisch erzeugten build-objekte.
* /results/ enthält alle von uns produzierten Graphen, Tabellen etc.
* /data/ enthält alle input und output files, die wichtig für den Encryption und Decryption Modus sind. 

## 3) Umsetzung des HPP-Automaten
Die Implementation aller HPP-Regeln befindet sich in /hpp_rules.hpp und /hpp_rules.cpp.
Das Grid ist standardmäßig als Torus angelegt (Randspalten und Randzeilen sind in direkter Nachbarschaft zueinander).

## 4) Parallelisierung
Wir verwenden MPI zur Zerlegung des Grids in 2D-Blöcke (domain-decomposition).
Aufgrund der Torus-Topologie des Grids besitzt jeder Block 2 Nachbarn (up/down und right/left).
Für solche Aufgaben stellt MPI die MPI_Cart Befehlsgruppe zur Verfügung.
Für den Übergang von "Teilchen" zwischen Domänen müssen diese miteinader kommunizieren. 
Dies geschieht über das Einführen von zusätzlichen Ghost-Zellen (auch Halo-Zellen) an den Rändern jedes 2D-Blocks.

Mit der Nutzung von MPI_Isend und MPI_Irecv kann die Kommunikation der Halo-Zellen überlappend mit der Berechnung der inneren (von den Halo-Zellen unabhängigen) 





Fußnote: MPI wurde auch zur Parallelsierung für einlesen und schreiben der Daten verwendet (siehe io.cpp), jedoch ist bei standard Gridgrößen und Nachrichtenlängen die Laufzeit zu > 98% im Haupt-Loop. 
Das Hauptaugemnerk liegt also auf der Parallelisierung dieses Loops.

## 5) Desing-Choices
- Wir speichern zusätzliche Metadaten (Originalgröße der Nachricht, Beginn der Nachricht im Grid, Gridgröße) in einem separaten File. Denkbar wäre auch die Verwendung eins "Stopp-bytes" um Ende (und ggf. Beginn) der Nachricht zu markieren.
- Die Gridgröße richtet sich standarmäßig nach der Größe der Nachricht. Feste Gridgrößen sind jedoch ebenfalls einstellbar und wurden zur Erzeugung der ERgebnisse verwendet.

## 6) Allgemeine Hinweise zum Code
-Das etwas längliche umrechnen und mappen (...)






