FC = gfortran
FFLAGS = -O3 -fopenmp -std=f2008 -fdefault-real-8 -fdefault-double-8

OBJS = \
./aerosol_mod.o \
./boundary_conditions_mod.o \
./chemistry_mod.o \
./derivatives_mod.o \
./drops_mod.o \
./dynamics_mod.o \
./grid_mod.o \
./main.o \
./meteorology_mod.o \
./parameterizations_mod.o \
./parameters_mod.o \
./prognostics_mod.o \
./radiation_mod.o \
./time_mod.o 

all: atmosphericModel

atmosphericModel: $(OBJS)
	@echo 'Building target: $@'
	$(FC) -fopenmp -o "atmosphericModel" $(OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o:  %.f90
	@echo 'Building file: $<'
	$(FC) -c $(FFLAGS) -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o:  %.f
	@echo 'Building file: $<'
	$(FC) -c $(FFLAGS) -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

aerosol_mod.o:  aerosol_mod.f90 parameters_mod.o time_mod.o

boundary_conditions_mod.o:  boundary_conditions_mod.f90 parameters_mod.o prognostics_mod.o

chemistry_mod.o:  chemistry_mod.f90 grid_mod.o parameters_mod.o time_mod.o

derivatives_mod.o:  derivatives_mod.f90 grid_mod.o

drops_mod.o:  drops_mod.f90 parameters_mod.o time_mod.o

dynamics_mod.o:  dynamics_mod.f90 aerosol_mod.o derivatives_mod.o grid_mod.o parameters_mod.o prognostics_mod.o time_mod.o


grid_mod.o:  grid_mod.f90 parameters_mod.o

main.o:  main.f90 aerosol_mod.o boundary_conditions_mod.o drops_mod.o dynamics_mod.o grid_mod.o parameterizations_mod.o parameters_mod.o prognostics_mod.o radiation_mod.o time_mod.o

meteorology_mod.o:  meteorology_mod.f90 grid_mod.o parameters_mod.o time_mod.o


parameterizations_mod.o:  parameterizations_mod.f90 aerosol_mod.o chemistry_mod.o parameters_mod.o prognostics_mod.o radiation_mod.o time_mod.o

parameters_mod.o:  parameters_mod.f90

prognostics_mod.o:  prognostics_mod.f90 aerosol_mod.o chemistry_mod.o derivatives_mod.o drops_mod.o grid_mod.o parameters_mod.o radiation_mod.o time_mod.o

radiation_mod.o:  radiation_mod.f90 parameters_mod.o time_mod.o

time_mod.o:  time_mod.f90 parameters_mod.o

clean:
	rm -rf $(OBJS) atmosphericModel
	-@echo ' '

.PHONY: all clean dependents


