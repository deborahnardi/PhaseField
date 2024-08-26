#pragma once

#include "Line.h"

class LineLoop
{
private:
    std::vector<Line *> lines;
    std::string entityName;
    int index;

public:
    LineLoop();
    LineLoop(const std::vector<Line *> &_lines, const int &_index = -1);
    ~LineLoop();

    // Setters
    void setLines(const std::vector<Line *> &_lines) { lines = _lines; }
    void setLine(const int &_index, Line *_line) { lines[_index] = _line; }
    void setIndex(const int &_index) { index = _index; }
    void setEntityName(const std::string &_entityName) { entityName = _entityName; }

    // Getters
    int getIndex() const { return index; }
    int getNumLines() const { return lines.size(); }
    std::string getEntityName() const { return entityName; }
    std::vector<Line *> setLines() { return lines; }
    Line *getLine(const int &_index) const { return lines[_index]; }
};