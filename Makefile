# Variables
CXX = g++
CXXFLAGS = -Wall -std=c++11
TARGET = main.cpp

# Default rule
all: small

# Compile the target
compile:
	$(CXX) $(CXXFLAGS) -o main.out $(TARGET)

# Clean rule
clean:
	rm -f *.out

# Run the program with parameters
facebook: compile
	./main.out facebook_combined.txt 4039 88234

random: compile
	./main.out random_1k_5k.txt 1000 9940

test: compile
	./main.out test.txt 11 25

.PHONY: all clean run