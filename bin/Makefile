F90 = gfortran
OBJ = lookup_table_module.o popiiistar_module.o testmodule.o
LIBS = 
FFLAGS = -frecord-marker=4 -O1 -ffree-line-length-none -x f95-cpp-input
VPATH = ../

testtable: $(OBJ)
	$(F90) $(OBJ) -o testmodule $(LIBS)
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
clean:	
	rm *.mod *.o
