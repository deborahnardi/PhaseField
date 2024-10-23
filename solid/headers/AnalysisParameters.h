#pragma once

#include <iostream>

#include "Point.h"
#include "Line.h"

class AnalysisParameters
{
private:
    int nSteps = 1., maxItNR = 1000.;
    double tol = 1.e-6;
    double ubar, deltaTime;
    double residStaggered = 1.e30;
    int maxItStaggered = 1000;
    std::vector<double> dispByStep;

public:
    AnalysisParameters();
    ~AnalysisParameters();

    int getNSteps() const { return nSteps; }
    int getMaxNewtonRaphsonIt() const { return maxItNR; }
    double getDeltaTime() const { return deltaTime; }
    double getTol() const { return tol; }
    double getResidStaggered() const { return residStaggered; }
    int getMaxItStaggered() const { return maxItStaggered; }

    void setNSteps(const int &_nSteps) { nSteps = _nSteps; }
    void setMaxNewtonRaphsonIt(const int &_maxItNR) { maxItNR = _maxItNR; }
    void setTol(const double &_tol) { tol = _tol; }
    void setResidStaggered(const double &_residStaggered) { residStaggered = _residStaggered; }
    void setMaxItStaggered(const int &_maxItStaggered) { maxItStaggered = _maxItStaggered; }
    void setDeltaTime(const double &_deltaTime) { deltaTime = _deltaTime; }
};