#pragma once

#include <iostream>

#include "Point.h"
#include "Line.h"

class AnalysisParameters
{
private:
    int nSteps, maxItNum = 1000;
    double tol = 1.e-8;
    double ubar;
    std::vector<double> dispByStep;

public:
    AnalysisParameters();
    ~AnalysisParameters();

    int getNSteps() const { return nSteps; }
    int getMaxItNum() const { return maxItNum; }
    double getTol() const { return tol; }

    void setNSteps(const int &_nSteps) { nSteps = _nSteps; }
    void setMaxItNum(const int &_maxItNum) { maxItNum = _maxItNum; }
    void setTol(const double &_tol) { tol = _tol; }
};