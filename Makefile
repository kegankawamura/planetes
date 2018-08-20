all:

test:

integrator_test: integrator_test.o NumIntegrators.o
	g++ -o bin/integrator_test obj/NumIntegrators.o obj/integrator_test.o
	
integrator_test.o: test/integrator_test.cpp
	g++ -o obj/integrator_test.o -c test/integrator_test.cpp

NumIntegrators.o: tools/NumIntegrators.cpp tools/NumIntegrators.h
	g++ -o obj/NumIntegrators.o -c tools/NumIntegrators.cpp

