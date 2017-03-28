################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../boundary_conditions_mod.f90 \
../chemistry_mod.f90 \
../derivatives_mod.f90 \
../dynamics_mod.f90 \
../grid_mod.f90 \
../main.f90 \
../meteorology_mod.f90 \
../parameterizations_mod.f90 \
../parameters_mod.f90 \
../prognostics_mod.f90 \
../radiation_mod.f90 \
../time_mod.f90 

OBJS += \
./boundary_conditions_mod.o \
./chemistry_mod.o \
./derivatives_mod.o \
./dynamics_mod.o \
./grid_mod.o \
./main.o \
./meteorology_mod.o \
./parameterizations_mod.o \
./parameters_mod.o \
./prognostics_mod.o \
./radiation_mod.o \
./time_mod.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

boundary_conditions_mod.o: ../boundary_conditions_mod.f90 parameters_mod.o prognostics_mod.o

chemistry_mod.o: ../chemistry_mod.f90 grid_mod.o parameters_mod.o time_mod.o

derivatives_mod.o: ../derivatives_mod.f90 grid_mod.o

dynamics_mod.o: ../dynamics_mod.f90 derivatives_mod.o grid_mod.o parameters_mod.o prognostics_mod.o time_mod.o

grid_mod.o: ../grid_mod.f90 parameters_mod.o

main.o: ../main.f90 boundary_conditions_mod.o dynamics_mod.o grid_mod.o parameterizations_mod.o parameters_mod.o prognostics_mod.o time_mod.o

meteorology_mod.o: ../meteorology_mod.f90 grid_mod.o parameters_mod.o time_mod.o

parameterizations_mod.o: ../parameterizations_mod.f90 chemistry_mod.o prognostics_mod.o radiation_mod.o time_mod.o

parameters_mod.o: ../parameters_mod.f90

prognostics_mod.o: ../prognostics_mod.f90 chemistry_mod.o derivatives_mod.o grid_mod.o parameters_mod.o radiation_mod.o time_mod.o

radiation_mod.o: ../radiation_mod.f90 parameters_mod.o time_mod.o

time_mod.o: ../time_mod.f90 parameters_mod.o


