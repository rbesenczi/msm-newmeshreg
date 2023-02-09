#include "meshregException.h"

namespace newmeshreg {

const char* MeshregException::what() const noexcept {
    std::cout << errmesg << std::endl;
    return errmesg;
}

} //namespace newmeshreg
