/*
            ======================================== MODE II CRACK PROPAGATION ============================================
            DECEMBER, 2024
            NO ENERGY SPLIT IMPLEMENTED YET
            DAMAGE LINE IMPOSED AS BOUNDARY CONDITION
            THE EXAMPLE IS AVALIABLE IN Ferreira, Marengo and Perego 2024 (https://doi.org/10.1016/j.cma.2024.117328)
            ==============================================================================================================
*/
std::string projectName = "phaseField2D-04";
Geometry *geo1 = new Geometry(projectName);
FEM *analysis1 = new FEM(projectName);
bool visualizeMesh = false;

PetscPrintf(PETSC_COMM_WORLD, "Running %s example...\n", projectName.c_str());

AnalysisParameters *params = new AnalysisParameters();
std::vector<Point *> points;
std::vector<Point *> pEllipse1;
std::vector<Point *> pEllipse2;
std::vector<Line *> lines;
std::vector<LineLoop *> lineLoops;
std::vector<Ellipse *> ellipses;
std::vector<PlaneSurface *> planeSurfaces;
std::vector<BoundaryCondition *> boundaryConditions;
std::vector<Material *> materials;

double L = 1.0;
double dL = 2e-4;
double elSize = 0.1 * L;
double ubar = 1e-3;

int np = 20;

// ================================ MESH GENERATION INFORMATION ================================
geo1->setAlgorithm(DELAUNAY);
geo1->setDimension(3);

points.push_back(geo1->addPoint({-L / 2, -L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({L / 2, -L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({L / 2, L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({-L / 2, L / 2, 0.0}, elSize));
points.push_back(geo1->addPoint({-L / 2, dL, 0.0}, elSize));
points.push_back(geo1->addPoint({0.0, 0.0, 0.0}, elSize));
points.push_back(geo1->addPoint({-L / 2, -dL, 0.0}, elSize));

double Lf = 0.4;
double alpha = 0.5;
double theta = 62. * M_PI / 180.;
double R1 = 0.5 * (1 + alpha) * Lf; // Major axis
double dc = Lf - R1;                // Distance from the center to the origin
double xc = dc * cos(theta);        // Center of the ellipse
double yc = dc * sin(theta);        // Center of the ellipse
double facR2 = 0.125;
double R2 = facR2 * R1; // Minor axis

// First point a = R1; b = R2
double A = pow(R2 * (cos(theta) - 2 * dL / L * sin(theta)), 2) + pow(R1 * (sin(theta) + 2 * dL / L * cos(theta)), 2);
double B = -2 * R2 * R2 * dc * (cos(theta) - 2 * dL / L * sin(theta));
double C = pow(R2 * dc, 2) - R1 * R1 * R2 * R2;

double x1 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
double y1 = -2 * dL / L * x1;

// Second point
A = pow(R2 * cos(theta), 2) + pow(R1 * sin(theta), 2);
B = -2 * R2 * R2 * dc * cos(theta);
C = pow(R2 * dc, 2) - R1 * R1 * R2 * R2;

double x2 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
double y2 = 0.0;

pEllipse1.push_back(geo1->addPoint({x1, -y1, 0.0}, elSize));

for (double thetai = -M_PI; thetai <= M_PI; thetai += (2 * M_PI) / np)
{
    double coord[] = {R1 * cos(thetai), R2 * sin(thetai)};           // Coordinates of the ellipse
    double xi = coord[0] * cos(theta) + coord[1] * sin(theta) + xc;  //
    double yi = -coord[0] * sin(theta) + coord[1] * cos(theta) - yc; //
    double coord2[] = {xi, yi};                                      // Rotating and translating the ellipse

    if (yi < 0.0)
        pEllipse1.push_back(geo1->addPoint({coord2[0], coord2[1], 0.0}, elSize));
}

pEllipse1.push_back(geo1->addPoint({x2, y2, 0.0}, elSize));
pEllipse2.push_back(pEllipse1.back());

for (double thetai = -M_PI; thetai <= M_PI; thetai += (2 * M_PI) / np)
{
    double coord[] = {R1 * cos(thetai), R2 * sin(thetai)};          // Coordinates of the ellipse
    double xi = coord[0] * cos(theta) - coord[1] * sin(theta) + xc; //
    double yi = coord[0] * sin(theta) + coord[1] * cos(theta) + yc; //
    double coord2[] = {xi, yi};                                     // Rotating and translating the ellipse

    if (yi > 0.0)
        pEllipse2.push_back(geo1->addPoint({coord2[0], coord2[1], 0.0}, elSize));
}

pEllipse2.push_back(geo1->addPoint({x1, y1, 0.0}, elSize));

lines.push_back(geo1->addLine({points[0], points[1]})); // 0
lines.push_back(geo1->addLine({points[1], points[2]})); // 1
lines.push_back(geo1->addLine({points[2], points[3]})); // 2
lines.push_back(geo1->addLine({points[3], points[4]})); // 3

lines.push_back(geo1->addLine({points[4], pEllipse2.back()})); // 4
lines.push_back(geo1->addLine({pEllipse2.back(), points[5]})); // 5
lines.push_back(geo1->addLine({points[5], pEllipse1[0]}));     // 6
lines.push_back(geo1->addLine({pEllipse1[0], points[6]}));     // 7
lines.push_back(geo1->addLine({points[6], points[0]}));        // 8

lines.push_back(geo1->addSpline(pEllipse1));
lines.push_back(geo1->addSpline(pEllipse2));

lineLoops.push_back(geo1->addLineLoop({lines[0], lines[1], lines[2], lines[3], lines[4], lines[10], lines[9], lines[7], lines[8]}));

lineLoops.push_back(geo1->addLineLoop({lines[10], lines[9], lines[6], lines[5]}));

planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[0]}));
planeSurfaces.push_back(geo1->addPlaneSurface({lineLoops[1]}));

boundaryConditions.push_back(geo1->addBoundaryCondition(lines[0], DIRICHLET, {{X, 0.0}, {Y, 0.0}}));
boundaryConditions.push_back(geo1->addBoundaryCondition(lines[2], DIRICHLET, {{X, ubar}, {Y, 0.0}}));

materials.push_back(geo1->addMaterial(210000., 0.3));
materials[0]->setGriffithCriterion(2.7);
materials[0]->setL0(0.01); // Internal lenght of Phase Field model, in mm;

planeSurfaces[0]->setAttributes(materials[0], 1., SOLID_ELEMENT);
planeSurfaces[1]->setAttributes(materials[0], 1., SOLID_ELEMENT);

double meshMinSizeGlobal = 1.e-4, meshMaxSizeGlobal = 0.1, meshSizeFactorGlobal = 1.0;
double meshMinSize = 1.e-2, meshMaxSize = 0.05, meshDistMin = 0.0001, meshDistMax = 0.002;

geo1->setSurfaceRefinement({lineLoops[1]}, meshMinSize, meshMaxSize, meshDistMin, meshDistMax);

// geo1->setGlobalMeshSize(meshMinSizeGlobal, meshMaxSizeGlobal, meshSizeFactorGlobal);
// geo1->setRefiningFieldCurves({lines[9], lines[10]}, 1);
// geo1->setThresholdRefinement(meshMinSize, meshMaxSize, meshDistMin, meshDistMax, 1, 2);
// geo1->setBackgroundMesh({2}, 3);

geo1->GenerateMeshAPI(visualizeMesh);

// ================================ FEM INFORMATION ================================
params->setDeltaTime(1);
params->setNSteps(146);
params->setSolverType(EIterative);

// Generating the loading vector

analysis1->setLoadingVector3(ubar, 146);

auto boundaryFunction = [](const std::vector<double> &coord, const double &pseudoTime, DOF *dof, const std::vector<double> &load)
{
    if (dof->getDOFType() == X)
        if (coord[1] == 0.5)
        {
            double val = load[pseudoTime];
            dof->setValue(val);
            dof->setControlledDOF();
        }
};
analysis1->setBoundaryFunction(boundaryFunction);
analysis1->setPrescribedDamageField(false);
// //   ********************************** FEM INFORMATION **********************************
params->setSolverType(EIterative);
params->setTolStaggered(1.e-4);
params->calculateReactionForces(true);
analysis1->setAnalysisParameters(params);
analysis1->readGeometry(projectName + ".mir");
analysis1->setPrintMatrix(false);
// analysis1->solveFEMProblem();
analysis1->solvePhaseFieldProblem();