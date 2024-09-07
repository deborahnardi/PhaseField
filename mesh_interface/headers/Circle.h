#pragma once

#include "Point.h"
#include "Material.h"
#include "../../enumclass.hpp"

class Circle
{
private:
    int index;
    double area;
    ElementType elementType;
    std::vector<Point *> points;
    std::string entityName;
    Material *material = nullptr;

public:
    Circle();
    Circle(const std::vector<Point *> &_points, const int &_index);
    ~Circle();

    // Setters
    void setPoints(const std::vector<Point *> &_points) { points = _points; }
    void setIndex(const int &_index) { index = _index; }
    void setEntityName(const std::string &_entityName) { entityName = _entityName; }
    void setMaterial(Material *_material) { material = _material; }
    void setAttributes(Material *_material, const double _area, const ElementType _elementType);

    // Getters
    int getIndex() const { return index; }
    double getArea() const { return area; }
    ElementType getElementType() const { return elementType; }
    std::vector<Point *> getPoints() const { return points; }
    Point *getPoint(const int &_index) const { return points[_index]; }
    std::string getEntityName() const { return entityName; }
    Material *getMaterial() const { return material; }
};