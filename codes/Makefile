OBJECTS = nr.o cosmomodule.o mcmc.o test.o
COMPILER = mpif90
FLAGS = -O3 
.PHONY = clean

test.output: main.exe
	./main.exe > test.output

main.exe: $(OBJECTS)
	$(COMPILER) $(FLAGS) $(OBJECTS) -o main.exe

%.o: %.f90
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm *.o *.mod *.exe *.output
