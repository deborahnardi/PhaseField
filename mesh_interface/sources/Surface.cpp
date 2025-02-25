#include "../headers/Surface.h"

Surface::Surface() {}
Surface::Surface(std::vector<int> _wireTags, const int &_index)
    : wireTags(_wireTags), index(_index), entityName("s" + std::to_string(index + 1)) {}
Surface::~Surface() {}

void Surface::setAttributes(Material *_material, const double _thickness, const ElementType _elementType)
{
    material = _material;
    thickness = _thickness;
    elementType = _elementType;
}