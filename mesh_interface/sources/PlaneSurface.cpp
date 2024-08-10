#include "../headers/PlaneSurface.h"

PlaneSurface::PlaneSurface() {}
PlaneSurface::PlaneSurface(LineLoop *lineLoop, const int index)
    : Surface(lineLoop, index) {}
PlaneSurface::~PlaneSurface() {}