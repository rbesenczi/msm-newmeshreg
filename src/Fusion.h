#ifndef NEWMESHREG_FUSION_H
#define NEWMESHREG_FUSION_H

#ifdef HAS_HOCR
#include "ELC/ELC.h"
#endif

#include "DiscreteCostFunction.h"
#include "FastPD.h"

#define NUM_SWEEPS 2
#define MAX_IMPROVEMENTS 0

typedef	double REAL;
struct UnaryData	{ REAL buffer[2]; };
struct PairData		{ REAL buffer[4]; };
struct TripletData	{ REAL buffer[8]; };

namespace ELCReduce { template<typename T> class PBF; }

template<typename T> class QPBO;

namespace newmeshreg {

enum Reduction { ELC_HOCR, ELC_APPROX, HOCR };

class Fusion {
public:
    template<typename OPTIMIZER> static void reduce_and_convert(ELCReduce::PBF<REAL>&, OPTIMIZER&, Reduction);
    static double optimize(const std::shared_ptr<DiscreteModel>& energy, Reduction reductionMode, bool debug = false);
};

} //namespace newmeshreg

#endif //NEWMESHREG_FUSION_H
