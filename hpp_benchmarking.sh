#!/bin/bash

# Erst sauberen Build starten (nur einmal nötig)
make clean
make BUILD=par

# Ergebnisdatei vorbereiten
echo "Prozesse | Laufzeit (s) " > hpp_benchmarking_results.txt

# Liste der Prozesszahlen
for NP in 2 4 6 8
do
    echo "Ausführen mit $NP Prozessen:"
    OUTPUT=$(make run BUILD=par NPROCS=$NP)

    RUNTIME=$(echo "$OUTPUT" | grep "Total runtime:" | sed -E 's/.*Total runtime:[[:space:]]*([0-9.]+).*/\1/')
    
    echo "$NP | $RUNTIME" >> hpp_benchmarking_results.txt
done