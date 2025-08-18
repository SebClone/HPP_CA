(englisch Version below)

Dokumentation für das Projekt "HPP-Automaton" im Modul parallel computing.
Von Samuel Orth und Sebastian Roth.

TldR: 

Inhaltverzeichnis:
0) Was, wie, warum?
1) Nutzung
2) Projektstuktur
3) Umsetzung des HPP-Automaten
4) Desing-Choices
5) Parallelisierung
6) Allgemeine Hinweise zum Code


1) Nutzung
Wir nutzen ein Makefile und eine config (siehe /include/config.hpp).
Eintellbar in der Config:
Der Modus (Encryption/Decryption) lässt sich über eine CLI Flag auswählen.
Dafür: make run NP=4 APP_ARG="--decrypt" oder

2) Projektstuktur
/include/ enthält alle Header Files
/src/ enthält alle .cpp Quellfiles
/build/ enthält alle Objekte, die (...)
/results/ enthält alle von uns produzierten Graphen, Tabellen etc.
/data/ ist ein 


4) Umsetzung des HPP-Automaten
Alle Regeln zur Umsetzung der Verschlüsselung befinden sich in /hpp_rules.hpp und /hpp_rules.cpp.
Das Grid ist ein Torus

6) Parallelisierung
Wir verwenden MPI beim (...)
Diese Maßnahmen, obowhl hilfreich, minimieren die Laufzeit jedoch nur wenig.
Hauptaugenmerk leigt auf der Paralleleisierung des Haup-Loops (siehe...)
Dafür machen wir eine Blockweise Zerlegung von 

8) Desing-Choices

6) Allgemeine Hinweise zum Code
-Das etwas längliche umrechnen und mappen 






