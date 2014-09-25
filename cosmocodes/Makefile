OBJECTS = nr.o cosmomodule.o mcmc.o fisher.o test.o
COMPILER = gfortran
FLAGS = -O3 
FINAL = main.exe
.PHONY = clean



test.output: $(FINAL)
	time ./$(FINAL)  #> test.output

$(FINAL): $(OBJECTS)
	$(COMPILER) $(FLAGS) $(OBJECTS) -o $(FINAL)

%.o: %.f90
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm *.o *.mod *.output
