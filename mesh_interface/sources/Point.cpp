#include "../headers/Point.h"

Point::Point() {}
Point::Point(const std::vector<double> &_coordinates, const double &_lc, const int &_index)
    : coordinates(_coordinates), lc(_lc), index(_index), entityName("p" + std::to_string(index + 1)) {}
Point::~Point() {}