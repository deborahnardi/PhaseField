#include "../headers/Line.h"

Line::Line() {}
Line::Line(const std::vector<Point *> &_points, const int &_index)
    : points(_points), index(_index), entityName("l" + std::to_string(index + 1)) {}
Line::~Line() {}

void Line::setAttributes(Material *_material, const double _area, const ElementType _elementType)
{
    material = _material;
    area = _area;
    elementType = _elementType;
}