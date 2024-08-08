#pragma once

enum SolverType
{
    EMumps,
    ESuiteSparse,
    EIterative
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

enum ElementName
{
    VERTEX = 15,
    // LINE ELEMENTS
    L2 = 1,
    L3 = 8,
    L4 = 26,
    NURBS_CURVE = 101,
    // SURFACE ELEMENTS
    T3 = 2,
    T6 = 9,
    T10 = 21,
    Q4 = 3,
    Q9 = 10,
    Q16 = 36,
    NURBS_SURFACE = 102,
    // VOLUME ELEMENTS
    TET4 = 4,
    TET10 = 11,
    TET20 = 29,
    HEX8 = 5,
    HEX27 = 12,
    HEX64 = 92,
    NURBS_VOLUME = 103
};

enum DOFType
{
    X = 0,
    Y = 1,
    Z = 2,
    DISP_X = 0,
    DISP_Y = 1,
    DISP_Z = 2,
    ROT_X = 3,
    ROT_Y = 4,
    ROT_Z = 5,
    VEC_X = 3,
    VEC_Y = 4,
    VEC_Z = 5,
    RATE_THICKNESS = 6,
    ALL
};

enum PlaneAnalysis
{
    PLANE_STRESS = 0,
    PLANE_STRAIN = 1
};

enum ElementType
{
    TRUSS_ELEMENT = 0,
    PLANE_ELEMENT = 1,
    SHELL_ELEMENT = 2,
    SOLID_ELEMENT = 3,
    NONE = -1
};