CXX       = g++
CXXFLAGS  = -std=c++17 -O2

SRCS      = main.cpp hpp_encryptor.cpp
OBJS      = $(SRCS:.cpp=.o)
TARGET    = encryptor

.PHONY: all run clean

all: $(TARGET)

run: $(TARGET)
	@./$(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

%.o: %.cpp hpp_encryptor.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
