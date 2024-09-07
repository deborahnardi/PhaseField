#include "../headers/Line.h"

Line::Line() {}
Line::Line(const std::vector<Point *> &_points, const int &_index)
    : points(_points), index(_index) {}
Line::~Line() {}

void Line::setAttributes(Material *_material, const double _area, const ElementType _elementType)
{
    material = _material;
    area = _area;
    elementType = _elementType;
}