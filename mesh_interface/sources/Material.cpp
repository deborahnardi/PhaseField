#include "../headers/Material.h"

Material::Material() {};
Material::Material(int _index, double _poisson, double _young, PlaneAnalysis _planeAnalysis)
    : index(_index), poisson(_poisson), young(_young), planeAnalysis(_planeAnalysis)
{
    shearModulus = young / (2. * (1. + poisson));

    switch (planeAnalysis)
    {
    case PLANE_STRAIN:
        lameConstant = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
        break;
    case PLANE_STRESS:
        lameConstant = young * poisson / (1. - poisson * poisson);
        break;
    }
}
Material::~Material() {};