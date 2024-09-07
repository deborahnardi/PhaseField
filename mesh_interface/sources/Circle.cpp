#include "../headers/Circle.h"

Circle::Circle() {}
Circle::Circle(const std::vector<Point *> &_points, const int &_index)
    : points(_points), index(_index) {}
Circle::~Circle() {}

void Circle::setAttributes(Material *_material, const double _area, const ElementType _elementType)
{
    material = _material;
    area = _area;
    elementType = _elementType;
}