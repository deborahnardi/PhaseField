#include "../headers/Surface.h"

Surface::Surface() {}
Surface::Surface(LineLoop *_lineLoop, const int &_index)
    : lineLoop(_lineLoop), index(_index) {}
Surface::~Surface() {}

void Surface::setAttributes(Material *_material, const double _thickness, const ElementType _elementType)
{
    material = _material;
    thickness = _thickness;
    elementType = _elementType;
}