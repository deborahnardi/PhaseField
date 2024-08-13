#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

class Node
{
private:
    int index;
    std::vector<double> initialCoordinates;

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
};