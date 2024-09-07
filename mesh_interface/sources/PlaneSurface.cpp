#include "../headers/PlaneSurface.h"

PlaneSurface::PlaneSurface() {}
PlaneSurface::PlaneSurface(std::vector<int> _wireTags, const int index)
    : Surface(_wireTags, index) {}
PlaneSurface::~PlaneSurface() {}