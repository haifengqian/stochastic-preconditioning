

SOURCES= main.cpp

OBJECTS= $(SOURCES:.cpp=.o)

default: testing

testing: $(OBJECTS)
	 g++ -m64 -o testing $(OBJECTS) -L../../src -lhybrid_solver 

$(OBJECTS): $(SOURCES)

$(OBJECTS): %.o: %.cpp
	 g++ -m64 -O3 -Wall -I../../src -c $< -o $@

clean:
	rm *.o
	rm testing
