#pragma once

#include "../../enumclass.hpp"
#include <string>

class Material
{
private:
    int index;
    double poisson, young, shearModulus, lameConstant, griffithCriterion, l0 = 0.;
    PlaneAnalysis planeAnalysis;
    std::string constitutiveType = "ISOTROPIC";

public:
    Material();
    Material(int _index, double _poisson, double _young, PlaneAnalysis _planeAnalysis = PLANE_STRESS);
    ~Material();

    int getIndex() { return index; };
    std::string getName() { return "Material_" + std::to_string(index + 1); }
    std::string getConstitutiveType() { return constitutiveType; }
    double getPoisson() { return poisson; }
    double getYoungModulus() { return young; }
    double getShearModulus() { return shearModulus; }
    double getLameConstant() { return lameConstant; }
    double getL0() const { return l0; }
    double getGriffithCriterion() { return griffithCriterion; }
    PlaneAnalysis getPlaneAnalysis() { return planeAnalysis; }

    void setIndex(int _index) { index = _index; }
    void setPoisson(double _poisson) { poisson = _poisson; }
    void setYoungModulus(double _young) { young = _young; }
    void setShearModulus(double _shearModulus) { shearModulus = _shearModulus; }
    void setPlaneAnalysis(PlaneAnalysis _planeAnalysis) { planeAnalysis = _planeAnalysis; }
    void setLameConstant(double _lameConstant) { lameConstant = _lameConstant; }
    void setGriffithCriterion(double _griff) { griffithCriterion = _griff; }
    void setL0(double _l0) { l0 = _l0; }
    void setConstitutiveType(std::string _constitutiveType) { constitutiveType = _constitutiveType; }

    void Lame(const double E[2][2], double S[2][2]);
};