std::string projectName = "square";
Geometry *geo1 = new Geometry(projectName, false); // true if has inclusions
bool visualizeMesh = false;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<LineLoop *> lineLoops;
std::vector<PlaneSurface *> planeSurfaces;

geo1->setAlgorithm(DELAUNAY);
geo1->setMeshSizeFactor(1.0);
geo1->setDimension(3);

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}, 0.25));
points.push_back(geo1->addPoint({1.0, 0.0, 0.0}, 0.25));
points.push_back(geo1->addPoint({1.0, 1.0, 0.0}, 0.25));
points.push_back(geo1->addPoint({0.0, 1.0, 0.0}, 0.25));

lines.push_back(geo1->addLine({points[0], points[1]}));
lines.push_back(geo1->addLine({points[1], points[2]}));
lines.push_back(geo1->addLine({points[2], points[3]}));
lines.push_back(geo1->addLine({points[3], points[0]}));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3]}));
planeSurfaces.push_back(geo1->addPlaneSurface(lineLoops[0]));

geo1->InitializeGmshAPI(visualizeMesh);