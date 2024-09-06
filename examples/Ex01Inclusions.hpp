std::string projectName = "rock_mesh1";
Geometry *geo1 = new Geometry(projectName, true); // true if has inclusions
FEM *solid1 = new FEM(projectName);
bool visualizeMesh = true;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

std::vector<Inclusion *> inclusions;
std::vector<MeshFactor *> factors;
std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<LineLoop *> lineLoops;
std::vector<PlaneSurface *> planeSurfaces;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

// geo1->setEdgeLength(1000.);
geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(2);
geo1->setMeshSizeFactor(1.0);

materials.push_back(geo1->addMaterial(1000., 0.2));
materials.push_back(geo1->addMaterial(2000., 0.2, INCLUSIONS));

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}));
points.push_back(geo1->addPoint({geo1->getEdgeLength(), 0.0, 0.0}));
points.push_back(geo1->addPoint({geo1->getEdgeLength(), geo1->getEdgeLength(), 0.0}));
points.push_back(geo1->addPoint({0.0, geo1->getEdgeLength(), 0.0}));

lines.push_back(geo1->addLine({points[0], points[1]}));
lines.push_back(geo1->addLine({points[1], points[2]}));
lines.push_back(geo1->addLine({points[2], points[3]}));
lines.push_back(geo1->addLine({points[3], points[0]}));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3]}));
planeSurfaces.push_back(geo1->addPlaneSurface(lineLoops[0]));

inclusions.push_back(geo1->addInclusion(0.05, 0.20, 30., 0.25, 0.25, 0.));
inclusions.push_back(geo1->addInclusion(0.05, 0.20, 60., 0.75, 0.75, 0.));
inclusions.push_back(geo1->addInclusion(0.10, 0.10, 90., 0.50, 0.50, 0.));
inclusions.push_back(geo1->addInclusion(0.10, 0.10, 45., 0.25, 0.75, 0.));

factors.push_back(geo1->addMeshFactor(0.1, 1.0, 1.2, 1e-4, 1e-1));

boundaryConditions.push_back(geo1->addBoundaryCondition(lines[3], DIRICHLET, {{X, 0.}, {Y, 0.}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[1], NEUMANN, {{X, 10.}}));

geo1->InitializeGmshAPI(visualizeMesh);

// solid1->readGeometry(projectName + ".mir");