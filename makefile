#-----------------------------------------------------#
# SIMMSUS - Makefile com perfis gfortran e ifx
# Objetos em build/obj e módulos em build/mod
#-----------------------------------------------------#

# Alvo padrão (pode ser sobrescrito via "make gfortran" ou "make ifx")
FC      ?= gfortran
FFLAGS  ?= -O2 -std=f2008 -fopenmp -fmax-stack-var-size=1073741824

EXEC = simmsus.ex

# Pastas de build
OBJDIR := build/obj
MODDIR := build/mod

# Fontes e objetos
SRC = variables.f90 subroutines.f90 input.f90 statistics.f90 simmsus.f90 main.f90
OBJ = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))

.PHONY: all clean gfortran ifx dirs

# Perfis de compilação
GFORTRAN_FLAGS = -O2 -std=f2008 -fopenmp -ffree-form -ffree-line-length-none -fmax-stack-var-size=1073741824 -J$(MODDIR) -I$(MODDIR)
IFX_FLAGS      = -O2 -qopenmp -heap-arrays -module $(MODDIR) -I$(MODDIR)

gfortran:
	$(MAKE) all FC=gfortran FFLAGS='$(GFORTRAN_FLAGS)'

ifx:
	$(MAKE) all FC=ifx FFLAGS='$(IFX_FLAGS)'

# Alvo principal
all: dirs $(EXEC)

# Binário
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

# Garante criação das pastas de build
dirs: $(OBJDIR) $(MODDIR)

$(OBJDIR) $(MODDIR):
	mkdir -p $@

# Regra genérica de compilação (gera .o em OBJDIR e .mod em MODDIR)
$(OBJDIR)/%.o: %.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

# Dependências entre módulos (garantem ordem correta)
$(OBJDIR)/subroutines.o: $(OBJDIR)/variables.o
$(OBJDIR)/input.o:       $(OBJDIR)/variables.o $(OBJDIR)/subroutines.o
$(OBJDIR)/statistics.o:  $(OBJDIR)/variables.o $(OBJDIR)/subroutines.o
$(OBJDIR)/simmsus.o:     $(OBJDIR)/variables.o $(OBJDIR)/subroutines.o
$(OBJDIR)/main.o:        $(OBJDIR)/variables.o $(OBJDIR)/subroutines.o $(OBJDIR)/input.o $(OBJDIR)/statistics.o $(OBJDIR)/simmsus.o

# Limpeza
clean:
	rm -rf $(OBJDIR) $(MODDIR) $(EXEC)

