#include "../headers/FEM.h"

void FEM::setReversibleDisp()
{
    /*
    In setReversibleDisp, the total applied displacement is divided into nSteps steps.
    The displacement is applied until the total displacement is reached.
    After that, the displacement is reversed until the initial position is reached.
    Then, the displacement is applied until the total displacement is reached again, but in the opposite direction.
    An unloading is perfomed again, until the initial position (in the opposite direction) is reached.
    */
    /*
        4 loading-unloading cycles
    */

    double ubar;

    for (auto dof : globalDOFs)
        if (dof->isControlledDOF())
        {
            ubar = dof->getDirichletValue(); // Considering that ONLY ONE DOF is controlled in the entire analysis
            break;
        }

    int nSteps = params->getNSteps();

    double disp = ubar / nSteps;

    for (int i = 0; i < nSteps; i++)
        dispByStep.push_back(disp * (i + 1));

    for (int i = 0; i < nSteps; i++)
        dispByStep.push_back((ubar + disp) - disp * (i + 1)); // Unloading

    for (int i = 0; i < nSteps; i++)
        dispByStep.push_back(-disp * (i + 1)); // Loading in the opposite direction

    for (int i = 0; i < nSteps; i++)
        dispByStep.push_back(-(ubar + disp) + disp * (i + 1)); // Loading in the opposite direction
}
void FEM::solvePhaseFieldProblem()
{
    DOF *controlledDOF;
    for (auto dof : globalDOFs)
        if (dof->isControlledDOF())
            controlledDOF = dof;

    for (int iStep = 0; iStep < dispByStep.size(); iStep++)
    {
        controlledDOF->setDirichletValue(dispByStep[iStep]);
    }
}