##################################################################

# simple user-application Makefile using the version of POOMA
# installed in the 'POOMADIR' directory.
##################################################################

### include the POOMA makefile stub, to get compiler flags and libraries

POOMADIR	=	/home/cli/lib/pooma241
POOMAARCH	=	linux

include $(POOMADIR)/$(POOMAARCH)/lib/Makefile.pooma


### the name of the example code to compile

MAIN 	= main
SUPPORT = Boundary.h utils.h functions.h Stencil.h PDESolver.h\
		  BurgerEquation.h ShallowWater.h LevelSet.h MultiVar.h\
		  LaxFriedrich.h Godunov.h EnoWeno.h
INCLUDE = -L/usr/lib64 -lnetcdf_c++
EXE		= run


### the main target for this makefile

# should include --openmp
$(EXE): $(MAIN).cpp $(SUPPORT)
	$(POOMA_CXX) -O3 --openmp -o $(EXE) $(MAIN).cpp $(POOMA_INCLUDES) $(POOMA_DEFINES) $(POOMA_LIBS) $(INCLUDE)

### clean things up a bit

clean:
	rm -rf $(MAIN) $(MAIN).o $(MAIN).ii ii_files ti_files $(EXE)
