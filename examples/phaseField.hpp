std::string projectName = "phaseField1D";
Geometry *geo1 = new Geometry(projectName); // true if has inclusions
FEM *analysis1 = new FEM(projectName);
bool visualizeMesh = true;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

AnalysisParameters *param = new AnalysisParameters();
std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

double L = 1.0;
double userNodes = 51;

geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

for (int i = 0; i < userNodes; i++)
    points.push_back(geo1->addPoint({L * i / (userNodes - 1), 0.0, 0.0}));

for (int i = 0; i < userNodes - 1; i++)
    lines.push_back(geo1->addLine({points[i], points[i + 1]}));

for (int i = 0; i < userNodes; i++)
    boundaryConditions.push_back(geo1->addBoundaryCondition(points[i], DIRICHLET, {{Y, 0.0}}));

boundaryConditions.push_back(geo1->addBoundaryCondition(points[0], DIRICHLET, {{X, 0.0}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(points[userNodes - 1], DIRICHLET, {{X, 0.005}}, true)); // true for controlled DOF

materials.push_back(geo1->addMaterial(21000., 0.));
materials[0]->setGriffithCriterion(0.35);
materials[0]->setL0(0.1); // Internal lenght of Phase Field model, in mm;

for (int i = 0; i < userNodes - 1; i++)
    lines[i]->setAttributes(materials[0], 1.0, TRUSS_ELEMENT);

materials.push_back(geo1->addMaterial(21000., 0.));
materials[1]->setGriffithCriterion(0.35 * 0.5);
materials[1]->setL0(0.1);

int flaw1 = ((userNodes - 1) / 2) - 1;
int flaw2 = (userNodes - 1) / 2;

lines[flaw1]->setAttributes(materials[1], 1.0, TRUSS_ELEMENT);
lines[flaw2]->setAttributes(materials[1], 1.0, TRUSS_ELEMENT);

param->setNSteps(200);

// ********************************** MESH GENERATION INFORMATION **********************************
geo1->GenerateMeshAPI(visualizeMesh);
// ********************************** FEM INFORMATION **********************************************
analysis1->readGeometry(projectName + ".mir");
analysis1->setAnalysisParameters(param);
// analysis1->setPrintMatrix(true);
// analysis1->solveFEMProblem();
//  analysis1->solveFEMProblemNoPetsc();
// ********************************** PHASE FIELD INFORMATION **************************************
analysis1->setReversibleDisp();
analysis1->solvePhaseFieldProblem();
