#pragma once

#include "Line.h"
#include <string>

class Wire
{
protected:
    int index;
    std::string entityName;

public:
    Wire();
    Wire(const int &_index);
    ~Wire();

    // Setters
    void setIndex(const int &_index) { index = _index; }
    void setEntityName(const std::string &_entityName) { entityName = _entityName; }

    // Getters
    int getIndex() const { return index; }
    std::string getEntityName() const { return entityName; }
};

class LineLoop : public Wire
{
private:
    std::vector<Line *> lines;

public:
    LineLoop();
    LineLoop(const std::vector<Line *> &_lines, const int &_index = -1);
    ~LineLoop();

    // Setters
    void setLines(const std::vector<Line *> &_lines) { lines = _lines; }
    void setLine(const int &_index, Line *_line) { lines[_index] = _line; }

    // Getters
    int getNumLines() const { return lines.size(); }
    std::vector<Line *> setLines() { return lines; }
    Line *getLine(const int &_index) const { return lines[_index]; }
};

class Ellipse : public Wire
{
private:
    double area, r1, r2, angle1, angle2;
    std::vector<double> center;
    std::vector<double> xAxis;

public:
    Ellipse();
    Ellipse(const std::vector<double> _points, double _r1, double _r2, const int _index, double _angle, double _angle2, std::vector<double> _xAxis);
    ~Ellipse();

    // Getters
    double getArea() const { return area; }
    std::vector<double> getCenter() const { return center; }
    double getR1() const { return r1; }
    double getR2() const { return r2; }
    double getAngle1() const { return angle1; }
    double getAngle2() const { return angle2; }
    std::vector<double> getXAxis() const { return xAxis; }
};

class EllipseArc : public Wire
{
private:
    int startTag, centerTag, majorTag, endTag;

public:
    EllipseArc();
    EllipseArc(const int _startTag, const int _centerTag, const int _majorTag, const int _endTag, const int _index);
    ~EllipseArc();

    // Getters
    double getStartTag() const { return startTag; }
    double getCenterTag() const { return centerTag; }
    double getMajorTag() const { return majorTag; }
    double getEndTag() const { return endTag; }
};