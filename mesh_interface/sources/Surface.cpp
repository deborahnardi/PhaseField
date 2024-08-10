#include "../headers/Surface.h"

Surface::Surface() {}
Surface::Surface(LineLoop* _lineLoop, const int &_index)
    : lineLoop(_lineLoop), index(_index) {}
Surface::~Surface() {}