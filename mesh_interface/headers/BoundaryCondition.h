#pragma once

#include <iostream>
#include <vector>
#include <string> // Add this line
#include "../../enumclass.hpp"

class BoundaryCondition
{
private:
    int index;
    bool controledDOF;
    std::string entityName;
    std::vector<std::pair<DOFType, double>> dofValues;
    BoundaryType bType;

public:
    BoundaryCondition();
    BoundaryCondition(int const &_index, const std::string &_entityName, const std::vector<std::pair<DOFType, double>> &_dofValues, const BoundaryType &_bType, const bool &_controledDOF);
    ~BoundaryCondition();

    int getIndex() const { return index; }
    bool isControledDOF() const { return controledDOF; }
    std::string getEntityName() const { return entityName; }
    std::vector<std::pair<DOFType, double>> getDOFValues() const { return dofValues; }
    BoundaryType getBType() const { return bType; }

    void setIndex(const int &_index) { index = _index; }
    void setEntityName(const std::string &_entityName) { entityName = _entityName; }
    void setDOFValues(const std::vector<std::pair<DOFType, double>> &_dofValues) { dofValues = _dofValues; }
    void setBType(const BoundaryType &_bType) { bType = _bType; }
};