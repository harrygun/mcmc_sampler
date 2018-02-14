CC= icc
FC= ifort



# C main() call fortran subroutine, for ifort.
COMFFLAGS = -nofor_main    
FLIB = -L/usr/local/lib  -lm   


MYLIBPATH= ~/documents/mylib
GENLIB   = lmylib
FISHLIB  = lfish




CFLAGS =  -g -I/usr/local/include -I $(MYLIBPATH)  -static
GSLIB = -L/usr/local/lib -lgsl -lgslcblas -lm


MYLIB = -L $(MYLIBPATH) -$(GENLIB) -$(FISHLIB)
#MYLIB = -L $(MYLIBPATH) -$(GENLIB) 


LIB  = $(GSLIB) $(MYLIB)




.SUFFIXES : .o .c
.c.o :
	$(CC) $(CFLAGS) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<





# workspace path
WORKSPACE = ../workspace



default: 




opt: 
	$(CC) $^ -o $(WORKSPACE)/$@  $(LIB)




clean:
	rm -f *.o *.c~
