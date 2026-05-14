# Fortran SP_Ace — PROGRAM space
#
# Build: make -C Fortran-95_code
# Produces ./space (Linux/macOS) or space.exe (MinGW on Windows)
#
# Use gfortran explicitly: GNU make predefines FC=f77 for legacy Fortran,
# so FC ?= gfortran never applies. Sources are ordered for MODULE deps.

FC      := gfortran
FFLAGS  ?= -O2 -ffree-line-length-none
LDFLAGS ?=

# Dependency order: each file must compile after its MODULE dependencies.
SRC := num_type.f95 stats.f95 share.f95 data_lib.f95 space_pars.f95 error.f95 \
	read_sp_ll.f95 read_GCOG.f95 func_poly.f95 make_model.f95 fit_cont.f95 \
	utils.f95 minimize.f95 uncertains2.f95 write_res.f95 interfaces.f95 space.f95

TARGET := space

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET) $(TARGET).exe *.o *.mod
