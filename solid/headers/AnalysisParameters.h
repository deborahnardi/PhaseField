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
    double residStaggered = 1.e30;
    int maxItStaggered = 1000;
    std::vector<double> dispByStep;

public:
    AnalysisParameters();
    ~AnalysisParameters();

    int getNSteps() const { return nSteps; }
    int getMaxItNum() const { return maxItNum; }
    double getTol() const { return tol; }
    double getResidStaggered() const { return residStaggered; }
    int getMaxItStaggered() const { return maxItStaggered; }

    void setNSteps(const int &_nSteps) { nSteps = _nSteps; }
    void setMaxItNum(const int &_maxItNum) { maxItNum = _maxItNum; }
    void setTol(const double &_tol) { tol = _tol; }
    void setResidStaggered(const double &_residStaggered) { residStaggered = _residStaggered; }
    void setMaxItStaggered(const int &_maxItStaggered) { maxItStaggered = _maxItStaggered; }
};