#include "../headers/Wire.h"

Wire::Wire() {}
Wire::Wire(const int &_index)
    : index(_index), entityName("w" + std::to_string(index + 1)) {}
Wire::~Wire() {}

LineLoop::LineLoop() {}
LineLoop::LineLoop(const std::vector<Line *> &_lines, const int &_index)
    : Wire(_index), lines(_lines) {}
LineLoop::~LineLoop() {}

Ellipse::Ellipse() {}
Ellipse::Ellipse(const std::vector<double> _points, double _r1, double _r2, const int _index, double _angle, double _angle2, std::vector<double> _xAxis)
    : Wire(_index), r1(_r1), r2(_r2), angle1(_angle), angle2(_angle2), center(_points), xAxis(_xAxis) {}
Ellipse::~Ellipse() {}