CC = g++
ODIR = obj


all:

clean:
	rm bin/*
	rm obj/*

test:

integrator_test: integrator_test.o NumIntegrators.o
	$(CC) -o bin/integrator_test obj/NumIntegrators.o obj/integrator_test.o
	
integrator_test.o: test/integrator_test.cpp
	$(CC) -o obj/integrator_test.o -c test/integrator_test.cpp

NumIntegrators.o: tools/NumIntegrators.cpp tools/NumIntegrators.h
	$(CC) -o obj/NumIntegrators.o -c tools/NumIntegrators.cpp

