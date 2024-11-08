#pragma once

#include <iostream>

#include "Point.h"
#include "Line.h"

class AnalysisParameters
{
private:
    int nSteps = 1, maxItNR = 1000, maxIterPSOR = 1000, maxItStaggered = 1000;
    double tolNR = 1.e-6, tolPSOR = 1.e-6, tolStaggered = 1.e-5;
    double ubar, deltaTime;
    double residStaggered = 1.e30, residPSOR = 1e20;
    std::vector<double> dispByStep;

public:
    AnalysisParameters();
    ~AnalysisParameters();

    int getNSteps() const { return nSteps; }
    int getMaxNewtonRaphsonIt() const { return maxItNR; }
    double getDeltaTime() const { return deltaTime; }
    double getTolNR() const { return tolNR; }
    double getResidStaggered() const { return residStaggered; }
    int getMaxItStaggered() const { return maxItStaggered; }
    int getResPSOR() const { return residPSOR; }
    int getMaxItPSOR() const { return maxIterPSOR; }
    double getTolPSOR() const { return tolPSOR; }
    double getTolStaggered() const { return tolStaggered; }

    void setNSteps(const int &_nSteps) { nSteps = _nSteps; }
    void setMaxNewtonRaphsonIt(const int &_maxItNR) { maxItNR = _maxItNR; }
    void setTolNR(const double &_tol) { tolNR = _tol; }
    void setResidStaggered(const double &_residStaggered) { residStaggered = _residStaggered; }
    void setMaxItStaggered(const int &_maxItStaggered) { maxItStaggered = _maxItStaggered; }
    void setDeltaTime(const double &_deltaTime) { deltaTime = _deltaTime; }
};