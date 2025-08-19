(englisch Version below)

#Dokumentation für das Projekt "HPP-Automaton" im Modul parallel computing.
Von Samuel Orth und Sebastian Roth.

##TldR: 

##Inhaltverzeichnis:
0. Projektinfo
1. Nutzung
2.  Projektstuktur
3.  Umsetzung des HPP-Automaten
4.  Desing-Choices
5.  Parallelisierung
6.  Allgemeine Hinweise zum Code

##0)
Dieses repository implementiert einen Lattice Gas Cellular Automaton - ein symmetrischer (...) nach dem Vorbild von (...).
Im Rahmen dieses Projekt wurde die encryption/decryption mithilfe von MPI und OpenMP parallelisisert

##1) Nutzung
Wir nutzen ein Makefile (/Makefile) und eine Konfigurationsdatei (/include/config.hpp, /src/config.cpp).
Einstellbar in der Config:
-iterations             (typ: integer, default: 1000, Aufgabe: )
-grid_size              (typ: integer, default: automatic, Aufgabe: )
-wall_density           (typ: double , default:, Aufgabe: )
-seed                   (typ: uint64_t, default: random, Aufgabe: )
-dump_frames            (typ: integer, default:, Aufgabe: )
-frame_interval         (typ: integer, default:, Aufgabe: )
-input                  (typ: integer, default:, Aufgabe: )
-enc_bin                (typ: integer, default:, Aufgabe: )
-meta                   (typ: integer, default:, Aufgabe: )
-key                    (typ: integer, default:, Aufgabe: )
-output                 (typ: integer, default:, Aufgabe: )

Einstellungen im Makefile:
-mpic++ Compiler
-build target heißt "hpp_mpi_app"
-basis flags für (...) -std=c++20 -Wall -Wextra -Wpedantic
-OpenMP:
-Laufzeitparameter: 
-CLI Flag:



Der Modus (Encryption/Decryption) lässt sich über eine CLI Flag auswählen.
Dafür: make run NP=4 APP_ARG="--decrypt" oder

##2) Projektstuktur
*/include/ enthält alle Header Files
*/src/ enthält alle .cpp Quellfiles
*/.build/ enthält alle automatisch erzeugten build-objekte
*/results/ enthält alle von uns produzierten Graphen, Tabellen etc.
*/data/ ist ein 


##3)
Die Implementation aller HPP-Regeln befindet sich in /hpp_rules.hpp und /hpp_rules.cpp.
Das Grid ist ein Torus

7) Parallelisierung
Wir verwenden MPI beim (...)
Diese Maßnahmen, obowhl hilfreich, minimieren die Laufzeit jedoch nur wenig.
Hauptaugenmerk leigt auf der Paralleleisierung des Haup-Loops (siehe...)
Dafür machen wir eine Blockweise Zerlegung von 

8) Desing-Choices

9) Allgemeine Hinweise zum Code
-Das etwas längliche umrechnen und mappen 






