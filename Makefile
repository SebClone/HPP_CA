# --- OpenMP automatic use for macOS---
# For macOS users: find latest g++ version and use it as OMP compiler (istall via homebrew if needed)
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

# Quellcode: alle .cpp im Verzeichnis (alte app_config.cpp ggf. ausschließen)
EXCLUDE_SRCS ?=
SRCS       := $(filter-out $(EXCLUDE_SRCS),$(wildcard *.cpp))
SRCS       := $(filter-out $(EXCLUDE_SRCS),$(wildcard src/*.cpp))
OBJS       := $(SRCS:.cpp=.o)
DEPS       := $(OBJS:.o=.d)

# Includes (anpassen/erweitern, falls du z. B. ein include/-Verzeichnis nutzt)
# CPPFLAGS   ?= -I.
CPPFLAGS   ?= -I. -Isrc -Iinclude

# C++-Standard und Warnungen
BASEFLAGS  := -std=c++20 -Wall -Wextra -Wpedantic

# Build-Typ: release|debug (make BUILD=debug)
BUILD      ?= release
ifeq ($(BUILD),debug)
  CXXFLAGS ?= -O0 -g -fno-omit-frame-pointer
  # Sanitizer optional (kommentiere aus, wenn nicht vorhanden/gewünscht)
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
  # Wenn OpenMP aus ist, unterdrücke ggf. Pragma-Warnungen:
  BASEFLAGS += -Wno-unknown-pragmas
endif

# Alles zusammenführen
CXXFLAGS  += $(BASEFLAGS) $(OMPFLAGS)
LDFLAGS   += $(OMPFLAGS)

# Optional: Architektur-Tuning (bei Bedarf)
# CXXFLAGS += -march=native


.PHONY: all clean run info

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Abhängigkeits-Dateien erzeugen (-MMD -MP)
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPS)

clean:
	$(RM) $(OBJS) $(DEPS) $(TARGET)

# Schnelltest: mpirun mit N Prozessen (make run NP=4)
NP ?= 4
run: $(TARGET)
	mpirun -np $(NP) ./$(TARGET)

# Zeigt aktuelle Einstellungen
info:
	@echo "CXX      = $(CXX)"
	@echo "OMPI_CXX = $(OMPI_CXX)"
	@echo "TARGET   = $(TARGET)"
	@echo "SRCS     = $(SRCS)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "LDFLAGS  = $(LDFLAGS)"
	@echo "CPPFLAGS = $(CPPFLAGS)"
	@echo "BUILD    = $(BUILD)"
	@echo "OMP      = $(OMP)"