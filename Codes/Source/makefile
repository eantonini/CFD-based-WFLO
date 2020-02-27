CC = g++

OBJS = PowerCalculation.o GradientCalculation.o CDGradientCalculation.o ReadInputs.o Geometry.o DiffusionProperties.o Simulation.o BCvelocity.o BCpressure.o BCk.o BCepsilon.o BComega.o BCnut.o ReadResults.o

all: $(OBJS)

PowerCalculation.o: PowerCalculation.cpp
	$(CC) -c PowerCalculation.cpp

GradientCalculation.o: GradientCalculation.cpp
	$(CC) -c GradientCalculation.cpp

CDGradientCalculation.o: CDGradientCalculation.cpp
	$(CC) -c CDGradientCalculation.cpp

ReadInputs.o: ReadInputs.cpp ReadInputs.h
	$(CC) -c ReadInputs.cpp

Geometry.o: Geometry.cpp Geometry.h
	$(CC) -c Geometry.cpp

DiffusionProperties.o: DiffusionProperties.cpp DiffusionProperties.h
	$(CC) -c DiffusionProperties.cpp

Simulation.o: Simulation.cpp Simulation.h
	$(CC) -c Simulation.cpp

BCvelocity.o: BCvelocity.cpp BCvelocity.h
	$(CC) -c BCvelocity.cpp

BCpressure.o: BCpressure.cpp BCpressure.h
	$(CC) -c BCpressure.cpp

BCk.o: BCk.cpp BCk.h
	$(CC) -c BCk.cpp

BCepsilon.o: BCepsilon.cpp BCepsilon.h
	$(CC) -c BCepsilon.cpp

BComega.o: BComega.cpp BComega.h
	$(CC) -c BComega.cpp

BCnut.o: BCnut.cpp BCnut.h
	$(CC) -c BCnut.cpp

ReadResults.o: ReadResults.cpp ReadResults.h
	$(CC) -c ReadResults.cpp

clean:
	rm $(OBJS)
