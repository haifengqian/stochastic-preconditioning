

SOURCES= HybridSolver.cpp \
         hsUtility.cpp \
         hsSymRmatCore.cpp \
         hsGenericRW.cpp

HEADERS= HybridSolver.h \
         hsUtility.h \
         hsSymRmatCore.h \
         hsGenericRW.h

OBJECTS= $(SOURCES:.cpp=.o)

default: libhybrid_solver.a

$(OBJECTS): $(SOURCES) $(HEADERS)

$(OBJECTS): %.o: %.cpp
	g++ -m64 -O3 -w -fpermissive -c $< -o $@

libhybrid_solver.a: $(OBJECTS)
	ar rcs libhybrid_solver.a $(OBJECTS)
	ranlib libhybrid_solver.a

clean:
	rm *.o
	rm *.a
