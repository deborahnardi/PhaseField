#pragma once

#include "../../enumclass.hpp"

class DOF
{
private:
    int index;
    DOFType type;
    double value;
    double vDirichlet = 0., vNeumann = 0.;
    bool bDirichlet = false, bNeumann = false;

public:
    DOF();
    DOF(const DOFType &_type, const double &_value = 0.);
    ~DOF();

    int getIndex() const { return index; }
    void setIndex(const int &_index) { index = _index; }

    DOFType getDOFType() const { return type; }
    void setDOFType(const DOFType &_type) { type = _type; }

    void setDirichlet() { bDirichlet = true; }
    void setDirichletValue(const double &_value) { vDirichlet = _value; }
    bool isDirichlet() const { return bDirichlet; }
    double getDirichletValue() const { return vDirichlet; }

    void setNeumann() { bNeumann = true; }
    void setNeumannValue(const double &_value) { vNeumann = _value; }
    bool isNeumann() const { return bNeumann; }
    double getNeumannValue() const { return vNeumann; }
};