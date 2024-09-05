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
    BOUNDARY_ELEMENT = 0,
    TRUSS_ELEMENT = 1,
    SOLID_ELEMENT = 2,
    NONE = -1
};

enum ApplyMaterial
{
    ALL = 0,
    INCLUSIONS = 1
};