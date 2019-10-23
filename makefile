CC = g++

EXE = PowerCalculation GradientCalculation CDGradientCalculation SLSQP

OBJS = Source/ReadInputs.o Source/Geometry.o Source/DiffusionProperties.o Source/Simulation.o Source/BCvelocity.o Source/BCpressure.o Source/BCk.o Source/BCepsilon.o Source/BComega.o Source/BCnut.o Source/ReadResults.o

all: $(EXE)

PowerCalculation:
	cd Source; make
	$(CC) -o PowerCalculation Source/PowerCalculation.o $(OBJS)

GradientCalculation:
	$(CC) -o GradientCalculation Source/GradientCalculation.o $(OBJS)

CDGradientCalculation:
	$(CC) -o CDGradientCalculation Source/CDGradientCalculation.o $(OBJS)

SLSQP: Source/SLSQP.cpp
	$(CC) -o SLSQP Source/SLSQP.cpp -lnlopt -lm -I/home/eantonini/NLopt/nlopt/include -L/home/eantonini/NLopt/nlopt/lib

clean:
	cd Source; make clean;
	rm $(EXE)
