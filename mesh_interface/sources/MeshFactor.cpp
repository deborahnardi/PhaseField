#include "../headers/Inclusion.h"

MeshFactor::MeshFactor(){}
MeshFactor::MeshFactor(const int &_index, const double &_meshMinFac, const double &_meshMaxFac, const double &_meshDistFac, const double &_meshMinSize, const double &_meshMaxSize)
    : index(_index), meshMinFac(_meshMinFac), meshMaxFac(_meshMaxFac), meshDistFac(_meshDistFac), meshMinSize(_meshMinSize), meshMaxSize(_meshMaxSize) {}
MeshFactor::~MeshFactor(){}