#pragma once

#include "../../enumclass.hpp"
#include <string>

class Material
{
private:
    int index;
    double poisson, young, shearModulus;
    std::string matType;
    ApplyMaterial whereToApply;
    PlaneAnalysis planeAnalysis;

public:
    Material();
    Material(int _index, double _poisson, double _young, PlaneAnalysis _planeAnalysis = PLANE_STRESS, const std::string &_matType = "ELASTIC", ApplyMaterial _whereToApply = ALL);
    ~Material();

    int getIndex() { return index; };
    std::string getName() { return "Material_" + std::to_string(index + 1); };
    std::string getMatType() { return matType; };
    double getPoisson() { return poisson; };
    double getYoungModulus() { return young; };
    double getShearModulus() { return shearModulus; };
    PlaneAnalysis getPlaneAnalysis() { return planeAnalysis; };

    void setIndex(int _index) { index = _index; };
    void setPoisson(double _poisson) { poisson = _poisson; };
    void setYoungModulus(double _young) { young = _young; };
    void setShearModulus(double _shearModulus) { shearModulus = _shearModulus; };
    void setName(const std::string &_name) { matType = _name; };
    void setMatType(const std::string &_matType) { matType = _matType; };
    void setPlaneAnalysis(PlaneAnalysis _planeAnalysis) { planeAnalysis = _planeAnalysis; };
};