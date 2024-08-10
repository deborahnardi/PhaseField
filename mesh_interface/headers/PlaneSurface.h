#include "Surface.h"

class PlaneSurface : public Surface
{
public:
    PlaneSurface();
    PlaneSurface(LineLoop *lineLoop, const int index);
    ~PlaneSurface();
};