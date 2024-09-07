std::string projectName = "squareEllipse";
Geometry *geo1 = new Geometry(projectName); // true if has inclusions
bool visualizeMesh = true;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<Ellipse *> ellipses;
std::vector<LineLoop *> lineLoops;
std::vector<PlaneSurface *> planeSurfaces;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

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

ellipses.push_back(geo1->addEllipse(0.5, 0.25, 30., {0.5, 0.5, 0.0}, 0.1));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[0]->getIndex(), ellipses[0]->getIndex()}));
planeSurfaces.push_back(geo1->addPlaneSurface({ellipses[0]->getIndex()}));

boundaryConditions.push_back(geo1->addBoundaryCondition(lines[3], DIRICHLET, {{X, 0.0}, {Y, 0.0}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[1], NEUMANN, {{X, 10.0}}));

materials.push_back(geo1->addMaterial(1000., 0.2));

planeSurfaces[0]->setAttributes(materials[0], 1., SOLID_ELEMENT);

geo1->InitializeGmshAPI(visualizeMesh);