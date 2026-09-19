CXX ?= c++
CXXFLAGS ?= -O3 -std=c++20
BOOST_ROOT ?=
INCLUDES := $(if $(BOOST_ROOT),-I$(BOOST_ROOT),)
TARGET := field_solver
SRC := main.cpp

all: $(TARGET)

$(TARGET): $(SRC) utils.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET) config.txt data/out.csv

clean:
	rm -f $(TARGET)

.PHONY: all run clean
