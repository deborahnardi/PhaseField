std::string projectName = "truss";
Geometry *geo1 = new Geometry(projectName);
FEM *truss = new FEM(projectName);
bool visualizeMesh = true;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<BoundaryCondition *> boundaryConditions;

double h = 1., b = 1., E = 1000., A0 = 1.;
geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}));
points.push_back(geo1->addPoint({b, h, 0.0}));
points.push_back(geo1->addPoint({2. * b, 0.0, 0.0}));

lines.push_back(geo1->addLine({points[0], points[1]}));
lines.push_back(geo1->addLine({points[1], points[2]}));

geo1->addTransfiniteLine({lines[0], lines[1]}, 1);

boundaryConditions.push_back(geo1->addBoundaryCondition(points[0], DIRICHLET, {{X, 0.}, {Y, 0.}, {Z, 0.}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(points[2], DIRICHLET, {{X, 0.}, {Y, 0.}, {Z, 0.}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(points[1], DIRICHLET, {{Z, 0.}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(points[1], NEUMANN, {{Y, -2000.}}));

geo1->GenerateMeshAPI(visualizeMesh);
// truss->readGeometry(projectName + ".mir");