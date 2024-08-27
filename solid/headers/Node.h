#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "../../enumclass.hpp"
#include "DOF.h"

class Node
{
private:
    int index;
    std::vector<double> initialCoordinates;
    std::vector<DOF *> dofs; // Degrees of freedom of the node

public:
    Node();
    Node(const int &_index, const std::vector<double> &_initialCoordinates);
    ~Node();

    // Setters
    void setIndex(const int &_index) { index = _index; }
    void setInitialCoordinates(const std::vector<double> &_initialCoordinates) { initialCoordinates = _initialCoordinates; }

    // Getters
    int getIndex() const { return index; }
    std::vector<double> getInitialCoordinates() const { return initialCoordinates; }
    double getX() const { return initialCoordinates[0]; }
    double getY() const { return initialCoordinates[1]; }
    double getZ() const { return initialCoordinates[2]; }

    void addDOF(DOF *_dof);
    std::vector<DOF *> getDOFs() const { return dofs; }
    DOF *getDOF(const int &_index) const { return dofs[_index]; }
    void setDOF(const int &_index, DOF *_dof) { dofs[_index] = _dof; }
    int getNumDOFs() const { return dofs.size(); }
};
