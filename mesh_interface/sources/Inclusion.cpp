#include "../headers/Inclusion.h"

Inclusion::Inclusion(){}
Inclusion::Inclusion(int _index, const double &_a, const double &_b, const double &_alpha, const double &_xc, const double &_yc)
    : index(_index), a(_a), b(_b), alpha(_alpha), xc(_xc), yc(_yc) {}
Inclusion::~Inclusion(){}