#pragma once

#include "Node.h"
#include "Material.h"

class BoundaryElement
{
public:
    BoundaryElement();
    BoundaryElement(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity);
    ~BoundaryElement();

    void getContribution() override {};
    void addCondition(BoundaryType _bdType, DOFType _type, double _value) override;
};