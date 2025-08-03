CXX       = mpic++ # Für MPI Einbindung
CXXFLAGS  = -std=c++17 -O2

SRCS      = main.cpp hpp_encryptor.cpp
OBJS      = $(SRCS:.cpp=.o)
TARGET    = encryptor

.PHONY: all run clean

all: $(TARGET)

run: $(TARGET)
	@mpirun -np 4 ./$(TARGET)   # hier mpirun zum Starten

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

%.o: %.cpp hpp_encryptor.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
