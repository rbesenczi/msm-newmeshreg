#ifndef NEWMESHREG_MESHREGEXCEPTION_H
#define NEWMESHREG_MESHREGEXCEPTION_H

#include <exception>
#include <iostream>

namespace newmeshreg {

class MeshregException : public std::exception {

public:
    const char* errmesg;

    explicit MeshregException(const char* msg) : errmesg(msg) {}

    const char* what() const noexcept override;
};

} //namespace newmeshreg

#endif //NEWMESHREG_MESHREGEXCEPTION_H
