#include "../headers/Line.h"

Line::Line() {}
Line::Line(const std::vector<Point *> &_points, const int &_index)
    : points(_points), index(_index) {}
Line::~Line() {}