std::string projectName = "phaseField2D";
Geometry *geo1 = new Geometry(projectName);
FEM *analysis1 = new FEM(projectName);
bool visualizeMesh = true;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

AnalysisParameters *params = new AnalysisParameters();
std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<LineLoop *> lineLoops;
std::vector<PlaneSurface *> planeSurfaces;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

double L = 1.0;
double dL = 0.0001;
double elSize = 0.1 * L;

// ================================ MESH GENERATION INFORMATION ================================
geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({L, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({L, L, 0.0}, elSize));
points.push_back(geo1->addPoint({0.0, L, 0.0}, elSize));
points.push_back(geo1->addPoint({0.0, L / 2 + dL, 0.0}, elSize));
points.push_back(geo1->addPoint({L / 2, L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({0.0, L / 2 - dL, 0.0}, elSize));

lines.push_back(geo1->addLine({points[0], points[1]}));
lines.push_back(geo1->addLine({points[1], points[2]}));
lines.push_back(geo1->addLine({points[2], points[3]}));
lines.push_back(geo1->addLine({points[3], points[4]}));
lines.push_back(geo1->addLine({points[4], points[5]}));
lines.push_back(geo1->addLine({points[5], points[6]}));
lines.push_back(geo1->addLine({points[6], points[0]}));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3], lines[4], lines[5], lines[6]}));

planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[0]}));

boundaryConditions.push_back(geo1->addBoundaryCondition(lines[0], DIRICHLET, {{X, 0.0}, {Y, 0.0}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[2], DIRICHLET, {{X, 0.0}, {Y, 0.004}}));

materials.push_back(geo1->addMaterial(210000., 0.3));
materials[0]->setGriffithCriterion(2.7);
materials[0]->setL0(0.01); // Internal lenght of Phase Field model, in mm;

planeSurfaces[0]->setAttributes(materials[0], 1., SOLID_ELEMENT);

double meshMinSizeGlobal = 1.e-4, meshMaxSizeGlobal = 0.1, meshSizeFactorGlobal = 1.0;
double meshMinSize = 0.005, meshMaxSize = 0.1, meshDistMin = 0.01, meshDistMax = 0.1;

geo1->setGlobalMeshSize(meshMinSizeGlobal, meshMaxSizeGlobal, meshSizeFactorGlobal);
geo1->setRefiningFieldCurves({lines[4], lines[5]}, 1); // lines 4 and 5 are the notch lines
geo1->setThresholdRefinement(meshMinSize, meshMaxSize, meshDistMin, meshDistMax, 1, 2);
geo1->setBoxRefinement(meshMinSize, meshMaxSize, L / 2, 1., 0.47, 0.53, 0.05, 3);
geo1->setBackgroundMesh({2, 3}, 4);

geo1->GenerateMeshAPI(visualizeMesh);

// ================================ FEM INFORMATION ================================
params->setDeltaTime(1);
params->setNSteps(80);
params->setSolverType(ESuiteSparse);

// Generating the loading vector

analysis1->setLoadingVector(0.004, 80);

auto boundaryFunction = [](const std::vector<double> &coord, const double &pseudoTime, DOF *dof, const std::vector<double> &load)
{
    if (dof->getDOFType() == Y)
        if (coord[1] == 1) // Applying load at the right end
        {
            double val = load[pseudoTime];
            dof->setValue(val);
            dof->setControlledDOF();
        }
};
analysis1->setBoundaryFunction(boundaryFunction);

// //   ********************************** FEM INFORMATION **********************************
params->setSolverType(EIterative);
params->setTolStaggered(1.e-4);
analysis1->setAnalysisParameters(params);
analysis1->readGeometry(projectName + ".mir");
analysis1->setPrintMatrix(false);
// analysis1->solveFEMProblem();
analysis1->solvePhaseFieldProblem();