#include "../headers/FEM.h"

void FEM::postProc()
{
    computeNodalStress();
}

void FEM::computeNodalStress()
{
    for (auto n : nodes)
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                n->setStress(i, j, 0.);
    for (auto el : elements)
        el->calculateStress();
}