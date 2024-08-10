/* 
This header file defines the Surface class, which is used to represent a surface in 1D, 2D or 3D space.
*/

#pragma once
#include "LineLoop.h"

class Surface
{
    protected: // not private so other subclasses can access it
        int index;
        std::string name;
        LineLoop* lineLoop;

    public:
        Surface();
        Surface(LineLoop* _lineLoop, const int &_index = -1);
        ~Surface();

        //Setters
        int getIndex() {return index;}
        LineLoop* getLineLoop() {return lineLoop;}
        
        // Getters
        void getIndex(int _index) {index = _index; }
        void setName(const std::string _name) {name = _name;}
        void setLineLoop(LineLoop *&_lineLoop) {lineLoop = _lineLoop;}
};