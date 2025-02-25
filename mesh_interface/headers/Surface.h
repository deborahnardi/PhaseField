/*
This header file defines the Surface class, which is used to represent a surface in 1D, 2D or 3D space.
*/

#pragma once
#include "Wire.h"
#include "Material.h"
#include "../../enumclass.hpp"

class Surface
{
protected: // not private so other subclasses can access it
    int index;
    double thickness;
    ElementType elementType;
    std::string name, entityName;
    std::vector<int> wireTags;
    Material *material = nullptr;

public:
    Surface();
    Surface(std::vector<int> _wireTags, const int &_index = -1);
    ~Surface();

    int getIndex() { return index; }
    double getThickness() { return thickness; }
    ElementType getElementType() { return elementType; }
    std::vector<int> getWireTag() { return wireTags; }
    std::string getName() { return name; }
    Material *getMaterial() { return material; }
    std::string getEntityName() { return entityName; }

    void setIndex(int _index) { index = _index; }
    void setName(const std::string _name) { name = _name; }
    void setEntityName(const std::string &_entityName) { entityName = _entityName; }
    void setMaterial(Material *_material) { material = _material; }
    void setAttributes(Material *_material, const double _thickness, const ElementType _elementType);
};