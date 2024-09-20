

SOURCES= generate.cpp \
         precondition.cpp

OBJECTS= $(SOURCES:.cpp=.o)

default: generate precondition

generate: $(OBJECTS)
	g++ -m64 -o generate generate.o -L../../src -lhybrid_solver 

precondition: $(OBJECTS)
	g++ -m64 -o precondition precondition.o -L../../src -lhybrid_solver 

$(OBJECTS): $(SOURCES)

$(OBJECTS): %.o: %.cpp
	g++ -m64 -O3 -w -I../../src -c $< -o $@

clean:
	rm *.o
	rm generate
	rm precondition

