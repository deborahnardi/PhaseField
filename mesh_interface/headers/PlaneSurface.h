#include "Surface.h"

class PlaneSurface : public Surface
{
private:
public:
    PlaneSurface();
    PlaneSurface(std::vector<int> _wireTags, const int index);
    ~PlaneSurface();

    void setEntityName(const std::string &_entityName) { entityName = _entityName; }
};