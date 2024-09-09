std::string projectName = "squareEllipse";
Geometry *geo1 = new Geometry(projectName); // true if has inclusions
FEM *analysis1 = new FEM(projectName);
bool visualizeMesh = true;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

std::vector<Point *> points;
std::vector<Line *> lines;
std::vector<Ellipse *> ellipses;
std::vector<LineLoop *> lineLoops;
std::vector<PlaneSurface *> planeSurfaces;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

double L = 1000.0;
double a0 = 0.05;

geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

points.push_back(geo1->addPoint({0.0, 0.0, 0.0}, 0.25 * L));
points.push_back(geo1->addPoint({L, 0.0, 0.0}, 0.25 * L));
points.push_back(geo1->addPoint({L, L, 0.0}, 0.25 * L));
points.push_back(geo1->addPoint({0.0, L, 0.0}, 0.25 * L));

lines.push_back(geo1->addLine({points[0], points[1]}));
lines.push_back(geo1->addLine({points[1], points[2]}));
lines.push_back(geo1->addLine({points[2], points[3]}));
lines.push_back(geo1->addLine({points[3], points[0]}));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3]}));

ellipses.push_back(geo1->addEllipse(a0 *L, 0.20, 30., {0.25 * L, 0.25 * L, 0.0 * L}));
ellipses.push_back(geo1->addEllipse(0.05 * L, 0.20, 60., {0.75 * L, 0.75 * L, 0.0 * L}));
ellipses.push_back(geo1->addEllipse(0.10 * L, 0.10, 90., {0.50 * L, 0.50 * L, 0.0 * L}));
ellipses.push_back(geo1->addEllipse(0.10 * L, 0.10, 45., {0.25 * L, 0.75 * L, 0.0 * L}));

planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[0], ellipses[0], ellipses[1], ellipses[2], ellipses[3]}));
planeSurfaces.push_back(geo1->addPlaneSurface({ellipses[0]}));
planeSurfaces.push_back(geo1->addPlaneSurface({ellipses[1]}));
planeSurfaces.push_back(geo1->addPlaneSurface({ellipses[2]}));
planeSurfaces.push_back(geo1->addPlaneSurface({ellipses[3]}));

boundaryConditions.push_back(geo1->addBoundaryCondition(lines[3], DIRICHLET, {{X, 0.0}, {Y, 0.0}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[1], NEUMANN, {{X, 10.0}}));

materials.push_back(geo1->addMaterial(1000., 0.2));
materials.push_back(geo1->addMaterial(5000., 0.15));

planeSurfaces[0]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[1]->setAttributes(materials[1], 1., SOLID_ELEMENT);
planeSurfaces[2]->setAttributes(materials[1], 1., SOLID_ELEMENT);
planeSurfaces[3]->setAttributes(materials[1], 1., SOLID_ELEMENT);
planeSurfaces[4]->setAttributes(materials[1], 1., SOLID_ELEMENT);

// ********************************** MESH GENERATION INFORMATION **********************************

double meshMinSizeIncl = 0.1 * a0 * L;
double meshMaxSizeIncl = 1.0 * a0 * L;
double meshDistMin = a0 * L;
double meshDistMax = 1.2 * a0 * L;

geo1->setGlobalMeshSize(1.e-4 * L, 1.e-1 * L, 10.);
geo1->setSurfaceRefinement({ellipses[0], ellipses[1], ellipses[2], ellipses[3]}, meshMinSizeIncl, meshMaxSizeIncl, meshDistMin, meshDistMax);
geo1->GenerateMeshAPI(visualizeMesh);

// ********************************** FEM INFORMATION **********************************
analysis1->readGeometry(projectName + ".mir");
