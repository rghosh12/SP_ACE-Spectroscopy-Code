# Fortran SP_Ace — PROGRAM space
# Build: make -C Fortran-95_code
# Produces ./space (Linux/macOS) or space.exe (MinGW on Windows)

FC      ?= gfortran
FFLAGS  ?= -O2 -ffree-line-length-none
LDFLAGS ?=

SRC := $(wildcard *.f95)
TARGET := space

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET) $(TARGET).exe
