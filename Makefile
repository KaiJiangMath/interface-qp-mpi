CXX 	= mpicxx
CFLAGS  = -DMPICH_IGNORE_CXX_SEEK -c -O2 -W -g
INCLUDE = -I./include -I/usr/include/suitesparse
LIBPATH = -L./lib
LDFLAGS = -lm -ldl -lfftw3_mpi -lfftw3 -lumfpack
OBJECTS	= MainCollect Data BasicOperators AuxInput FftwToolkit StableBulkPhase StableBulkParallel \
		  GenJacPoly ProjFourGJP CommonRotateProjBoxMat CommonFourGJP SysInitVal \
		  SparseOperators IterPrepare IterMethod DisplayResults MemFree

all:interface
# 建立动态链接
./lib/lib%.so:./src/%.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(INCLUDE) $(LDFLAGS)

interface:$(addsuffix .so, $(addprefix ./lib/lib, $(OBJECTS))) interface.cpp
		$(CXX) -g $(INCLUDE) interface.cpp -o $@ $(LIBPATH) $(addprefix -l, $(OBJECTS)) $(LDFLAGS)
clean: 
	-rm ./lib/*.so
	-rm interface
cleanout:
	-rm ./lib/*.so
	-rm interface
	-rm -r ./result
	-rm -r ./fig
	-rm -r ./vtk
	-rm -r ./plan
	-rm -r ./log
	-rm -r ./parallel
