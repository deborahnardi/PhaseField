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

void Material::Lame(const double E[2][2], double S[2][2])
{
    const double &trE = E[0][0] + E[1][1];
    S[0][0] = 2. * shearModulus * E[0][0] + lameConstant * trE;
    S[1][1] = 2. * shearModulus * E[1][1] + lameConstant * trE;
    S[0][1] = 2. * shearModulus * E[0][1];
    S[1][0] = 2. * shearModulus * E[1][0];
}