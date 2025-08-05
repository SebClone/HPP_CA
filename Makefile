NPROCS ?= $(shell getconf _NPROCESSORS_ONLN)

BUILD ?= seq

ifeq ($(BUILD),par)
  CXX     := mpicxx
  SRC     := main_parallel.cpp
  RUN_CMD := mpirun -np $(NPROCS) ./encrypt
else
  CXX     := g++
  SRC     := main.cpp
  RUN_CMD := ./encrypt
endif

CXXFLAGS := -O3 -std=c++17

.PHONY: all run clean
all: encrypt

encrypt:
	$(CXX) $(CXXFLAGS) $(SRC) -o encrypt

run: encrypt
	@echo "Running with BUILD=$(BUILD), NPROCS=$(NPROCS)"
	$(RUN_CMD)

clean:
	rm -f encrypt
