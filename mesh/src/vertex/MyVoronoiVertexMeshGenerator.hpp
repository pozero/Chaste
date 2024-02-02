#ifndef MYVORONOIVERTEXMESHGENERATOR_HPP_
#define MYVORONOIVERTEXMESHGENERATOR_HPP_

#include <vector>

#include "MutableVertexMesh.hpp"
#include "fortune_voronoi/FortuneAlgorithm.h"

class MyVoronoiVertexMeshGenerator
{
private:
    MutableVertexMesh<2, 2>* mpMesh;

public:
    MyVoronoiVertexMeshGenerator(std::vector<Vector2> const& sites, Box const& boudning_box);

    MyVoronoiVertexMeshGenerator() = default;

    ~MyVoronoiVertexMeshGenerator();

    MutableVertexMesh<2, 2>* GetMesh();
};

#endif
