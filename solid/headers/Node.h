#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include <petscsnes.h>
#include <petscksp.h>
#include <petscdraw.h>
#include <petscmat.h>
#include <metis.h>

#include <algorithm>

#include "../../enumclass.hpp"
#include "DOF.h"

class Node
{
private:
    int index, physicalEntity;
    std::vector<double> initialCoordinates;
    std::vector<DOF *> dofs; // Degrees of freedom of the node
    bool isDiscritized = false;

    std::vector<int> inverseIncidence;

    double stress[3][3] = {};

public:
    Node();
    Node(const int &_index, const std::vector<double> &_initialCoordinates);
    ~Node();

    // Setters
    void setIndex(const int &_index) { index = _index; }
    void setInitialCoordinates(const std::vector<double> &_initialCoordinates) { initialCoordinates = _initialCoordinates; }

    // Getters
    int getIndex() const { return index; }
    int getPhysicalEntity() const { return physicalEntity; }
    std::vector<double> getInitialCoordinates() const { return initialCoordinates; }
    double getX() const { return initialCoordinates[0]; }
    double getY() const { return initialCoordinates[1]; }
    double getZ() const { return initialCoordinates[2]; }
    double getIsDiscritized() const { return isDiscritized; }

    void addDOF(DOF *_dof);
    std::vector<DOF *> getDOFs() const { return dofs; }
    DOF *getDOF(const int &_index) const { return dofs[_index]; }
    void setDOF(const int &_index, DOF *_dof) { dofs[_index] = _dof; }
    int getNumDOFs() const { return dofs.size(); }
    void setIsDiscritized() { isDiscritized = true; }
    void setPhysicalEntity(const int &physicalEntity_) { physicalEntity = physicalEntity_; }

    void addInverseIncidence(int el);
    std::vector<int> getInverseIncidence() { return inverseIncidence; }

    void setStress(const int &i, const int &j, const double &val) { stress[i][j] = val; }
    void incrementStress(const int &i, const int &j, const double &val) { stress[i][j] += val; }
    double getStress(const int &i, const int &j) { return stress[i][j]; }
};
