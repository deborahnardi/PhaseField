#include "../headers/Material.h"

Material::Material() {};
Material::Material(int _index, double _poisson, double _young, PlaneAnalysis _planeAnalysis, const std::string &_matType, ApplyMaterial _whereToApply)
    : index(_index), poisson(_poisson), young(_young), planeAnalysis(_planeAnalysis), matType(_matType), whereToApply(_whereToApply)
{
    shearModulus = young / (2. * (1. + poisson));
}
Material::~Material() {};