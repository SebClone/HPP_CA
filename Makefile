# --- OpenMP automatic use for macOS---
# For macOS users: find latest g++ version and use it as OMP compiler (install via homebrew if needed)
ifeq ($(shell uname),Darwin)
  DETECTED_GCC := $(shell command -v g++-14 || command -v g++-13 || command -v g++-12)
  ifneq ($(DETECTED_GCC),)
    export OMPI_CXX=$(DETECTED_GCC)
  endif
endif

#
# --- Ensure CXX is MPI wrapper if default or c++/clang++ from environment ---
ifeq ($(origin CXX),default)
  CXX := mpic++
endif
ifeq ($(origin CXX),environment)
  ifeq ($(CXX),c++)
    CXX := mpic++
  endif
  ifeq ($(CXX),clang++)
    CXX := mpic++
  endif
endif

# -------------------------------------------
# --- Basis-Konfig ---
CXX        ?= mpic++
TARGET     ?= hpp_mpi_app

# Build-Ordner für .o/.d:
BUILD_DIR  ?= .build

# Quellcode: alle .cpp im Root und unter src/ (EXCLUDE_SRCS optional ausschließen)
EXCLUDE_SRCS ?=
SRCS := $(filter-out $(EXCLUDE_SRCS),$(wildcard *.cpp) $(wildcard src/*.cpp))

# Objekt-/Dep-Dateien im BUILD_DIR, Verzeichnisstruktur wird gespiegelt
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# Includes
CPPFLAGS   ?= -I. -Isrc -Iinclude
# --- Fix für nervige OpenMPI/MPICH C++-Bindings-Warnungen ---
# Deaktiviert die alten C++-Bindings (wir nutzen die C-API mit mpi.h).
CPPFLAGS   += -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX
# Falls du die C++-Bindings wirklich brauchst (MPI::), nimm die Zeile wieder raus.
# Alternativ (nicht so sauber): Warnung unterdrücken
# CXXFLAGS  += -Wno-cast-function-type

# C++-Standard und Warnungen
BASEFLAGS  := -std=c++20 -Wall -Wextra -Wpedantic

# Build-Typ: release|debug (make BUILD=debug)
BUILD      ?= release
ifeq ($(BUILD),debug)
  CXXFLAGS ?= -O0 -g -fno-omit-frame-pointer
  # Sanitizer optional:
  # CXXFLAGS += -fsanitize=address,undefined
  # LDFLAGS  += -fsanitize=address,undefined
else
  CXXFLAGS ?= -O3
endif

# OpenMP (make OMP=0 zum Deaktivieren)
OMP        ?= 1
ifeq ($(OMP),1)
  OMPFLAGS := -fopenmp
else
  OMPFLAGS :=
  BASEFLAGS += -Wno-unknown-pragmas
endif

# Alles zusammenführen
CXXFLAGS  += $(BASEFLAGS) $(OMPFLAGS)
LDFLAGS   += $(OMPFLAGS)

# Optional: Architektur-Tuning
# CXXFLAGS += -march=native

.PHONY: all clean run info rerun

all: $(TARGET)

# Link-Schritt: Binary bleibt im Projektroot
$(TARGET): $(OBJS)
	@mkdir -p $(@D)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Kompilieren: .o/.d ins BUILD_DIR (inkl. Unterordner), .d via -MF explizit setzen
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -MF $(@:.o=.d) -c $< -o $@

-include $(DEPS)

clean:
	$(RM) -r $(BUILD_DIR) $(TARGET)

# --- Laufzeitparameter ---
NP       ?= 4          # Anzahl MPI-Prozesse
RUN_ARGS ?=            # z.B. --bind-to core --map-by core
APP_ARGS ?=            # z.B. --decrypt oder -m encrypt

# Start: mpirun mit RUN_ARGS, Programm bekommt APP_ARGS
run: $(TARGET)
	mpirun -np $(NP) $(RUN_ARGS) ./$(TARGET) $(APP_ARGS)

# Komfort: clean + run
rerun: clean run

# Zeigt aktuelle Einstellungen
info:
	@echo "CXX       = $(CXX)"
	@echo "OMPI_CXX  = $(OMPI_CXX)"
	@echo "TARGET    = $(TARGET)"
	@echo "BUILD_DIR = $(BUILD_DIR)"
	@echo "SRCS      = $(SRCS)"
	@echo "OBJS      = $(OBJS)"
	@echo "CXXFLAGS  = $(CXXFLAGS)"
	@echo "LDFLAGS   = $(LDFLAGS)"
	@echo "CPPFLAGS  = $(CPPFLAGS)"
	@echo "BUILD     = $(BUILD)"
	@echo "OMP       = $(OMP)"
	@echo "NP        = $(NP)"
	@echo "RUN_ARGS  = $(RUN_ARGS)"
	@echo "APP_ARGS  = $(APP_ARGS)"