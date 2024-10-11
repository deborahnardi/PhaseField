#include "../headers/AnalysisParameters.h"

AnalysisParameters::AnalysisParameters() {}
AnalysisParameters::~AnalysisParameters() {}

/*
    In setReversibleDisp, the total applied displacement is divided into nSteps steps.
    The displacement is applied until the total displacement is reached.
    After that, the displacement is reversed until the initial position is reached.
    Then, the displacement is applied until the total displacement is reached again, but in the opposite direction.
    An unloading is perfomed again, until the initial position (in the opposite direction) is reached.
*/

void AnalysisParameters::setReversibleDisp(const double &_ubar, const Point *point)
{
    ubar = _ubar;

    double dispArray[nSteps * 4]{}; // 4 loading-unloading cycles
    double disp = ubar / nSteps;

    for (int i = 0; i < nSteps; i++)
        dispArray[i] = disp * (i + 1);

    for (int i = 0; i < nSteps; i++)
        dispArray[i + nSteps] = (ubar + disp) - disp * (i + 1); // Unloading

    for (int i = 0; i < nSteps; i++)
        dispArray[i + 2 * nSteps] = -disp * (i + 1); // Loading in the opposite direction

    for (int i = 0; i < nSteps; i++)
        dispArray[i + 3 * nSteps] = -(ubar + disp) + disp * (i + 1); // Loading in the opposite direction

    for (int i = 0; i < nSteps * 4; i++)
        dispByStep.push_back(dispArray[i]);
}