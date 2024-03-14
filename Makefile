all: manifold manifold-qgb paramtest

LIBGEOM=../libgeom
TRANSFINITE=../transfinite
QGB=../qgb

INCLUDES=-I$(LIBGEOM) -I$(TRANSFINITE)/src -I$(QGB)
LIBS=-lOpenMeshCore \
	-L$(TRANSFINITE)/release/transfinite -ltransfinite \
	-L$(LIBGEOM)/release -lgeom

CXXFLAGS=-std=c++20 -Wall -pedantic -DNDEBUG -O3 $(INCLUDES)

manifold: manifold.o param.o blend.o
	$(CXX) -o $@ $^ $(LIBS)

manifold-qgb: manifold-qgb.o param.o blend.o
	$(CXX) -o $@ $^ $(QGB)/qgb.o $(LIBS)

paramtest: paramtest.o
	$(CXX) -o $@ $< $(CXXFLAGS) $(LIBS)
