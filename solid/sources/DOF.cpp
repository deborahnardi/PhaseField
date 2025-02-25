#include "../headers/DOF.h"

DOF::DOF() {}
DOF::DOF(const DOFType &_type, const double &_value)
    : type(_type), value(_value) {}
DOF::~DOF() {}