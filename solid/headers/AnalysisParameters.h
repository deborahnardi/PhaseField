#pragma once

#include <iostream>

#include "Point.h"
#include "Line.h"

class AnalysisParameters
{
private:
    int nSteps = 1, maxItNR = 1000, maxIterPSOR = 1000, maxItStaggered = 100000;
    double tolNR = 1.e-6, tolPSOR = 1.e-6, tolStaggered = 1.e-5;
    double ubar, deltaTime;
    double residStaggered = 1.e30, residPSOR = 1e20;
    int maxIterEIterative = 100000;
    double tolEIterative = 1.e-8;
    std::vector<double> dispByStep;
    SolverType solverType = EMumps;
    bool calcReactionForces = false;
    std::string PFmodel = "AT2", reactionDir = "X";

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
    SolverType getSolverType() const { return solverType; }
    int getMaxIterEIterative() const { return maxIterEIterative; }
    double getTolEIterative() const { return tolEIterative; }
    bool getCalculateReactionForces() const { return calcReactionForces; }
    std::string getPFModel() const { return PFmodel; }
    std::string getReactionDir() const { return reactionDir; }

    void setNSteps(const int &_nSteps) { nSteps = _nSteps; }
    void setMaxNewtonRaphsonIt(const int &_maxItNR) { maxItNR = _maxItNR; }
    void setTolNR(const double &_tol) { tolNR = _tol; }
    void setResidStaggered(const double &_residStaggered) { residStaggered = _residStaggered; }
    void setMaxItStaggered(const int &_maxItStaggered) { maxItStaggered = _maxItStaggered; }
    void setDeltaTime(const double &_deltaTime) { deltaTime = _deltaTime; }
    void setSolverType(const SolverType &_solverType) { solverType = _solverType; }
    void setTolStaggered(const double &_tolStaggered) { tolStaggered = _tolStaggered; }
    void setCalculateReactionForces(const bool &_calcReactionForces) { calcReactionForces = _calcReactionForces; }
    void calculateReactionForces(const bool &_calcReactionForces) { calcReactionForces = _calcReactionForces; }
    void setReactionDir(const std::string &_reactionDir) { reactionDir = _reactionDir; }
    void setPFModel(const std::string &_PFmodel) { PFmodel = _PFmodel; }
};