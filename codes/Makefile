OBJECTS = cosmomodule.o test.o
COMPILER = gfortran
FLAGS = -O3 -fopenmp
.PHONY = clean

test.output: main.exe
	./main.exe > test.output

main.exe: $(OBJECTS)
	$(COMPILER) $(FLAGS) $(OBJECTS) -o main.exe

%.o: %.f90
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm *.o *.mod *.exe *.output
