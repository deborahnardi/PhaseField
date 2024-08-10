#pragma once

#include "Point.h"

class Line
{
private:
    std::vector<Point *> points;
    int index;

public:
    Line();
    Line(const std::vector<Point *> &_points, const int &_index);
    ~Line();

    //Setters
    void setPoints(const std::vector<Point *> &_points) { points = _points; }
    void setIndex(const int &_index) { index = _index; }

    //Getters
    int getIndex() const { return index; }

    std::vector<Point*> getPoints() const { return points; }
    Point *getPoint(const int &_index) const { return points[_index]; }
};
