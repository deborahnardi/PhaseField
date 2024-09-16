#pragma once

enum SolverType
{
    Sequential,
    Parallel
};

enum BoundaryType
{
    DIRICHLET = 0,
    NEUMANN = 1
};

enum MeshAlgorithm
{
    ADAPT = 1,
    AUTO = 2,
    DELAUNAY = 5,
    FRONT = 6,
    BAMG = 7,
    QUAD = 8,
    PACK = 9,
    HXT = 10
};

enum DOFType
{
    X = 0,
    Y = 1,
    Z = 2
};

enum PlaneAnalysis
{
    PLANE_STRESS = 0,
    PLANE_STRAIN = 1
};

enum ElementType
{
    TRUSS_ELEMENT = 0,
    SOLID_ELEMENT = 1,
    NONE = -1
};

enum DeformType
{
    DEFORM_NONE = 0,
    DEFORM_SHEAR = 1,
    DEFORM_STEP = 2,
    NUM_DEFORM_TYPES = 3
};

enum SolutionType
{
    SOL_VLAP_QUADRATIC = 0,
    SOL_ELAS_QUADRATIC = 1,
    SOL_VLAP_TRIG = 2,
    SOL_ELAS_TRIG = 3,
    SOL_ELAS_AXIAL_DISP = 4,
    SOL_ELAS_UNIFORM_STRAIN = 5,
    SOL_ELAS_GE = 6,
    SOL_MASS_QUADRATIC = 7,
    NUM_SOLUTION_TYPES = 8
};