#pragma once

#include "../../enumclass.hpp"

class DOF
{
private:
    int index = -1.;
    DOFType type;
    double value = 0.; // Value is the current value of the DOF
    double vDirichlet = 0., vNeumann = 0., vDamage = 0.;
    bool bDirichlet = false, bNeumann = false, isControlled = false, bDamage = false;

public:
    DOF();
    DOF(const DOFType &_type, const double &_value = 0.);
    ~DOF();

    int getIndex() const { return index; }
    double getValue() const { return value; }
    double getDamageValue() const { return vDamage; }
    void setIndex(const int &_index) { index = _index; }
    void setValue(const double &_value) { value = _value; }
    void incrementValue(const double &_increment) { value += _increment; }
    void setControlledDOF() { isControlled = true; }

    DOFType getDOFType() const { return type; }
    void setDOFType(const DOFType &_type) { type = _type; }

    void setDirichlet() { bDirichlet = true; }
    void setDirichletValue(const double &_value) { vDirichlet = _value, value = _value; }
    void setDamageValue(const double &_value) { vDamage = _value; }
    bool isDirichlet() const { return bDirichlet; }
    bool isControlledDOF() const { return isControlled; }
    double getDirichletValue() const { return vDirichlet; }

    void setNeumann() { bNeumann = true; }
    void setNeumannValue(const double &_value) { vNeumann = _value, value = _value; }
    bool isNeumann() const { return bNeumann; }
    double getNeumannValue() const { return vNeumann; }

    void setPrescribedDamage() { bDamage = true; }
};