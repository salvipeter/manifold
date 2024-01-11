all: manifold

FLAGS=-std=c++20 -Wall -pedantic -g -DDEBUG -O0
#FLAGS=-std=c++20 -Wall -pedantic -DNDEBUG -O3

LIBGEOM=../libgeom
TRANSFINITE=../transfinite

INCLUDES=-I$(LIBGEOM) -I$(TRANSFINITE)/src
LIBS=-lOpenMeshCore \
	-L$(TRANSFINITE)/release/transfinite -ltransfinite \
	-L$(LIBGEOM)/release -lgeom

manifold: manifold.cc
	$(CXX) -o $@ $< $(FLAGS) $(INCLUDES) $(LIBS)
