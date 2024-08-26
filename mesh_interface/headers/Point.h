#pragma once

#include <iostream>
#include <vector>
#include <string>

class Point
{
private:
    int index;
    std::vector<double> coordinates;
    std::string entityName;
    double lc; // characteristic length

public:
    Point();
    Point(const std::vector<double> &_coordinates, const double &_lc, const int &_index = -1.0);
    ~Point();

    // Setters
    void setCoordinates(const std::vector<double> &_coordinates) { coordinates = _coordinates; }
    void setLC(const double &_lc) { lc = _lc; }
    void setIndex(const int &_index) { index = _index; }
    void setEntityName(const std::string &_entityName) { entityName = _entityName; }

    // Getters
    int getIndex() const { return index; }
    std::vector<double> getCoordinates() const { return coordinates; }
    std::string getEntityName() const { return entityName; }
    double getX() const { return coordinates[0]; }
    double getY() const { return coordinates[1]; }
    double getZ() const { return coordinates[2]; }
    double getLC() const { return lc; }
};