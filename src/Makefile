include ${FSLCONFDIR}/default.mk

PROJNAME    = NewMeshReg
SOFILES     = libfsl-newmeshreg.so
LIBS        = -lfsl-newresampler -lfsl-newimage -lfsl-miscmaths -lfsl-utils
USRCXXFLAGS = -fopenmp
USRINCFLAGS = -DHAS_HOCR

OBJS       = reg_tools.o meshregException.o histogram2D.o \
             featurespace.o similarities.o DiscreteCostFunction.o \
             DiscreteModel.o Fusion.o  mesh_registration.o \
             rigid_costfunction.o FastPD.o graph.o

all: ${SOFILES}

libfsl-newmeshreg.so: ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}
