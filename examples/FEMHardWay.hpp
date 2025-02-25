/*
            ======================================== MODE I CRACK PROPAGATION ============================================
            FEB. 24th, 2025
            "HARD WAY" OF FEM MEANS THAT THERE ARE NO LOOPS OVER THE ELEMENTS ANYMORE, BUT OVER THE NODES
            ==============================================================================================================
*/
std::string projectName = "FEMHardWay";
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

double L = 1.0, E = 1000., A0 = 1.;
double elSize = 2 * L;
// ================================ MESH GENERATION INFORMATION ================================
geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({L, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({2 * L, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({3 * L, 0.0, 0.0}, elSize));

points.push_back(geo1->addPoint({0.0, L, 0.0}, elSize));
points.push_back(geo1->addPoint({L, L, 0.0}, elSize));
points.push_back(geo1->addPoint({2 * L, L, 0.0}, elSize));
points.push_back(geo1->addPoint({3 * L, L, 0.0}, elSize));

lines.push_back(geo1->addLine({points[0], points[1]})); // 0
lines.push_back(geo1->addLine({points[1], points[4]})); // 1
lines.push_back(geo1->addLine({points[4], points[0]})); // 2
lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[0]}));

lines.push_back(geo1->addLine({points[1], points[5]})); // 3
lines.push_back(geo1->addLine({points[5], points[4]})); // 4
lines.push_back(geo1->addLine({points[4], points[1]})); // 5
lineLoops.push_back(geo1->addLineLoop({lines[3], lines[4], lines[5]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[1]}));

lines.push_back(geo1->addLine({points[1], points[6]})); // 6
lines.push_back(geo1->addLine({points[6], points[5]})); // 7
lines.push_back(geo1->addLine({points[5], points[1]})); // 8
lineLoops.push_back(geo1->addLineLoop({lines[6], lines[7], lines[8]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[2]}));

lines.push_back(geo1->addLine({points[1], points[2]})); // 9
lines.push_back(geo1->addLine({points[2], points[6]})); // 10
lines.push_back(geo1->addLine({points[6], points[1]})); // 11
lineLoops.push_back(geo1->addLineLoop({lines[9], lines[10], lines[11]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[3]}));

lines.push_back(geo1->addLine({points[2], points[3]})); // 12
lines.push_back(geo1->addLine({points[3], points[6]})); // 13
lines.push_back(geo1->addLine({points[6], points[2]})); // 14
lineLoops.push_back(geo1->addLineLoop({lines[12], lines[13], lines[14]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[4]}));

lines.push_back(geo1->addLine({points[3], points[7]})); // 15
lines.push_back(geo1->addLine({points[7], points[6]})); // 16
lines.push_back(geo1->addLine({points[6], points[3]})); // 17
lineLoops.push_back(geo1->addLineLoop({lines[15], lines[16], lines[17]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[5]}));

boundaryConditions.push_back(geo1->addBoundaryCondition(points[0], DIRICHLET, {{X, 0.}, {Y, 0.}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(points[3], DIRICHLET, {{Y, 0.}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(points[7], NEUMANN, {{Y, -1.}}));

materials.push_back(geo1->addMaterial(E, 0.0));

planeSurfaces[0]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[1]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[2]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[3]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[4]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[5]->setAttributes(materials[0], 1., SOLID_ELEMENT);

geo1->GenerateMeshAPI(visualizeMesh);
analysis1->setAnalysisParameters(params);
analysis1->readGeometry(projectName + ".mir");
analysis1->setPrintMatrix(false);
analysis1->solveFEMProblem();