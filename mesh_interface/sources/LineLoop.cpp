#include "../headers/LineLoop.h"

LineLoop::LineLoop() {}
LineLoop::LineLoop(const std::vector<Line *> &_lines, const int &_index)
    : lines(_lines), index(_index) {}
LineLoop::~LineLoop() {}