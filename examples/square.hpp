std::string projectName = "square";
Geometry *geo1 = new Geometry(projectName); // true if has inclusions
FEM *analysis1 = new FEM(projectName);
bool visualizeMesh = false;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<LineLoop *> lineLoops;
std::vector<PlaneSurface *> planeSurfaces;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

double L = 1.0;

geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

double elSize = 0.5 * L;

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}, elSize)); // Those are lists (dynamic arrays) of pointers to the objects
points.push_back(geo1->addPoint({L, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({L, L, 0.0}, elSize));
points.push_back(geo1->addPoint({0.0, L, 0.0}, elSize));

lines.push_back(geo1->addLine({points[0], points[1]}));
lines.push_back(geo1->addLine({points[1], points[2]}));
lines.push_back(geo1->addLine({points[2], points[3]}));
lines.push_back(geo1->addLine({points[3], points[0]}));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3]}));

planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[0]}));

boundaryConditions.push_back(geo1->addBoundaryCondition(lines[3], DIRICHLET, {{X, 0.0}, {Y, 0.0}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[1], DIRICHLET, {{X, 1.0}}));

materials.push_back(geo1->addMaterial(1., 0.3));

planeSurfaces[0]->setAttributes(materials[0], 1., SOLID_ELEMENT);

// ********************************** MESH GENERATION INFORMATION **********************************

geo1->GenerateMeshAPI(visualizeMesh);

// ********************************** FEM INFORMATION **********************************
analysis1->readGeometry(projectName + ".mir");
analysis1->setPrintMatrix(true);
analysis1->solveFEMProblem();
// analysis1->solveFEMProblemNoPetsc();
//  ********************************* NUMERICAL INTEGRATION *********************************
