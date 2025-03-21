/*
            ======================================== MODE I CRACK PROPAGATION ============================================
            DECEMBER 5th, 2024
            DAMAGE LINE IMPOSED AS BOUNDARY CONDITION
            THE EXAMPLE IS AVALIABLE IN Ferreira, Marengo and Perego 2024 (https://doi.org/10.1016/j.cma.2024.117328)
            ==============================================================================================================
*/
std::string projectName = "phaseField2D-03";
Geometry *geo1 = new Geometry(projectName);
FEM *analysis1 = new FEM(projectName);
bool visualizeMesh = false;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

AnalysisParameters *params = new AnalysisParameters();
std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<LineLoop *> lineLoops;
std::vector<PlaneSurface *> planeSurfaces;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

double L = 1.0;
double elSize = 0.1 * L;
double ubar = 4e-4;

// ================================ MESH GENERATION INFORMATION ================================
geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({L, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({L, L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({L / 2, L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({0.0, L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({L, L, 0.0}, elSize));
points.push_back(geo1->addPoint({0.0, L, 0.0}, elSize));

lines.push_back(geo1->addLine({points[0], points[1]}));
lines.push_back(geo1->addLine({points[1], points[2]}));
lines.push_back(geo1->addLine({points[2], points[3]}));
lines.push_back(geo1->addLine({points[3], points[4]}));
lines.push_back(geo1->addLine({points[4], points[0]}));
lines.push_back(geo1->addLine({points[2], points[5]}));
lines.push_back(geo1->addLine({points[5], points[6]}));
lines.push_back(geo1->addLine({points[6], points[4]}));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3], lines[4]}));
lineLoops.push_back(geo1->addLineLoop({lines[5], lines[6], lines[7], lines[3], lines[2]}));

planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[0]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[1]}));

boundaryConditions.push_back(geo1->addBoundaryCondition(lines[0], DIRICHLET, {{X, 0.0}, {Y, 0.0}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[6], DIRICHLET, {{X, 0.0}, {Y, ubar}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[3], DAMAGE, {{D, 1.0}}));

geo1->setTractionBoundary({lines[6]});

materials.push_back(geo1->addMaterial(210000., 0.3));
materials[0]->setGriffithCriterion(2.7);
materials[0]->setL0(0.01); // Internal lenght of Phase Field model, in mm;

planeSurfaces[0]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[1]->setAttributes(materials[0], 1., SOLID_ELEMENT);

double meshMinSizeGlobal = 1, meshMaxSizeGlobal = 1, meshSizeFactorGlobal = 1;
// double meshMinSize = 0.05, meshMaxSize = 0.1, meshDistMin = 0.01, meshDistMax = 0.05;
double meshMinSize = 1, meshMaxSize = 1, meshDistMin = 1, meshDistMax = 1;

geo1->setGlobalMeshSize(meshMinSizeGlobal, meshMaxSizeGlobal, meshSizeFactorGlobal);
// geo1->setRefiningFieldCurves({lines[3]}, 1);
// geo1->setThresholdRefinement(meshMinSize, meshMaxSize, meshDistMin, meshDistMax, 1, 2);
// geo1->setBoxRefinement(meshMinSize, meshMaxSize, L / 2, 1., 0.48, 0.52, 0.05, 3);
// geo1->setBoxRefinement(meshMinSize, meshMaxSize, L / 2, 1., 0.47, 0.53, 0.05, 3);
geo1->setBackgroundMesh({2, 3}, 4);

geo1->GenerateMeshAPI(visualizeMesh);

// ================================ FEM INFORMATION ================================
params->setDeltaTime(1);
params->setNSteps(80);

// Generating the loading vector

analysis1->setLoadingVector2(ubar, 80);

auto boundaryFunction = [](const std::vector<double> &coord, const double &pseudoTime, DOF *dof, const std::vector<double> &load)
{
    if (dof->getDOFType() == Y)
        if (coord[1] == 1)
        {
            double val = load[pseudoTime];
            dof->setValue(val);
            dof->setControlledDOF();
        }
};
analysis1->setBoundaryFunction(boundaryFunction);
analysis1->setPrescribedDamageField(false);
//  //   ********************************** FEM INFORMATION **********************************
params->setSolverType(EMumps);
params->setTolStaggered(1.e-4);
params->calculateReactionForces(true);
params->setReactionDir("Y");
params->setPFModel("AT2");
analysis1->setAnalysisParameters(params);
analysis1->readGeometry(projectName + ".mir");
analysis1->setPrintMatrix(false);
// analysis1->solveFEMProblem();
analysis1->solvePhaseFieldProblem();