all: manifold manifold-normal manifold-qgb

#FLAGS=-std=c++20 -Wall -pedantic -g -DDEBUG -O0
FLAGS=-std=c++20 -Wall -pedantic -DNDEBUG -O3

LIBGEOM=../libgeom
TRANSFINITE=../transfinite
QGB=../qgb

INCLUDES=-I$(LIBGEOM) -I$(TRANSFINITE)/src -I$(QGB)
LIBS=-lOpenMeshCore \
	-L$(TRANSFINITE)/release/transfinite -ltransfinite \
	-L$(LIBGEOM)/release -lgeom

manifold: manifold.cc
	$(CXX) -o $@ $< $(FLAGS) $(INCLUDES) $(LIBS)

manifold-normal: manifold-normal.cc
	$(CXX) -o $@ $< $(FLAGS) $(INCLUDES) $(LIBS)

manifold-qgb: manifold-qgb.cc
	$(CXX) -o $@ $< $(QGB)/qgb.o $(FLAGS) $(INCLUDES) $(LIBS)
