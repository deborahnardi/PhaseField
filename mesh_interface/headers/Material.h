#pragma once

#include "../../enumclass.hpp"
#include <string>

class Material
{
private:
    int index;
    double poisson, young, shearModulus, lameConstant;
    PlaneAnalysis planeAnalysis;

public:
    Material();
    Material(int _index, double _poisson, double _young, PlaneAnalysis _planeAnalysis = PLANE_STRESS);
    ~Material();

    int getIndex() { return index; };
    std::string getName() { return "Material_" + std::to_string(index + 1); }
    double getPoisson() { return poisson; }
    double getYoungModulus() { return young; }
    double getShearModulus() { return shearModulus; }
    double getLameConstant() { return lameConstant; }
    PlaneAnalysis getPlaneAnalysis() { return planeAnalysis; }

    void setIndex(int _index) { index = _index; }
    void setPoisson(double _poisson) { poisson = _poisson; }
    void setYoungModulus(double _young) { young = _young; }
    void setShearModulus(double _shearModulus) { shearModulus = _shearModulus; }
    void setPlaneAnalysis(PlaneAnalysis _planeAnalysis) { planeAnalysis = _planeAnalysis; }
    void setLameConstant(double _lameConstant) { lameConstant = _lameConstant; }
};