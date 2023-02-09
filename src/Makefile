include ${FSLCONFDIR}/default.mk

PROJNAME = NewMeshReg
SOFILES  = libfsl-newmeshreg.so
LIBS     = -lfsl-fastPD2 -lfsl-newimage -lfsl-miscmaths -lfsl-utils -lfsl-newresampler
USRCXXFLAGS+= -DHAS_HOCR -fopenmp

all: ${SOFILES}

libfsl-newmeshreg.so: reg_tools.o meshregException.o histogram2D.o featurespace.o similarities.o DiscreteCostFunction.o DiscreteModel.o Fusion.o mesh_registration.o rigid_costfunction.o
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}
