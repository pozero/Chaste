#ifndef MYVORONOIVERTEXMEHSGENERATOR_HPP_
#define MYVORONOIVERTEXMEHSGENERATOR_HPP_

#include <cmath>
#include <vector>
#include <map>

#include "MutableVertexMesh.hpp"

class MyVoronoiVertexMeshGenerator
{
protected:
    MutableVertexMesh<2, 2>* mpMesh;

public:
    MyVoronoiVertexMeshGenerator(
        std::vector<double>& cell_target_areas);

    ~MyVoronoiVertexMeshGenerator();

    MutableVertexMesh<2, 2>* GetMesh();
};

#endif
