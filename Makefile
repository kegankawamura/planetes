CC = g++
ODIR = obj
SDIR = src
TESTDIR = $(SDIR)/test


all:

clean:
	rm bin/*
	rm $(ODIR)/*

test:

integrator_test: integrator_test.o NumIntegrators.o
	$(CC) -o bin/integrator_test $(ODIR)/NumIntegrators.o $(ODIR)/integrator_test.o
	
integrator_test.o: $(TESTDIR)/integrator_test.cpp
	$(CC) -o $(ODIR)/integrator_test.o -c $(TESTDIR)/integrator_test.cpp

NumIntegrators.o: src/tools/NumIntegrators.cpp src/tools/NumIntegrators.h
	$(CC) -o $(ODIR)/NumIntegrators.o -c src/tools/NumIntegrators.cpp

