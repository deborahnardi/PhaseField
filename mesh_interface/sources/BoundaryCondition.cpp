#include "../headers/BoundaryCondition.h"

BoundaryCondition::BoundaryCondition() {}
BoundaryCondition::BoundaryCondition(int const &_index, const std::string &_entityName, const std::vector<std::pair<DOFType, double>> &_dofValues, const BoundaryType &_bType)
    : index(_index), entityName(_entityName), dofValues(_dofValues), bType(_bType) {}
BoundaryCondition::~BoundaryCondition() {}