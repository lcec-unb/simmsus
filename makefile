#-----------------------------------------------------#
# SIMMSUS - Makefile com perfis gfortran e ifx
#-----------------------------------------------------#

# Alvo padrão (pode ser sobrescrito via "make gfortran" ou "make ifx")
FC      ?= gfortran
FFLAGS  ?= -O2 -std=f2008 -fopenmp -fmax-stack-var-size=1073741824

EXEC = simmsus.ex

SRC = variables.f90 subroutines.f90 input.f90 statistics.f90 simmsus.f90 main.f90
OBJ = $(SRC:.f90=.o)

.PHONY: all clean gfortran ifx

# Perfis de compilação
GFORTRAN_FLAGS = -O2 -std=f2008 -fopenmp -fmax-stack-var-size=1073741824
IFX_FLAGS      = -O2 -qopenmp -heap-arrays

gfortran:
	$(MAKE) all FC=gfortran FFLAGS='$(GFORTRAN_FLAGS)'

ifx:
	$(MAKE) all FC=ifx FFLAGS='$(IFX_FLAGS)'

# Alvo principal
all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

# Regra genérica de compilação
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Dependências entre módulos (garantem ordem correta)
subroutines.o: variables.o
input.o:       variables.o subroutines.o
statistics.o:  variables.o subroutines.o
simmsus.o:     variables.o subroutines.o
main.o:        variables.o subroutines.o input.o statistics.o simmsus.o

# Limpeza
clean:
	rm -f *.o *.mod $(EXEC)

