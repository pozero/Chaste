#include "MyVoronoiVertexMeshGenerator.hpp"

MyVoronoiVertexMeshGenerator::MyVoronoiVertexMeshGenerator(std::vector<Vector2> const& sites, Box const& boudning_box)
{
    FortuneAlgorithm fortune_algorithm{ sites };
    fortune_algorithm.construct();
    fortune_algorithm.bound(boudning_box);
    VoronoiDiagram voronoi_diagram = fortune_algorithm.getDiagram();

    unsigned vertex_idx = 0;
    std::map<const VoronoiDiagram::Vertex*, unsigned> vertex_indices{};
    std::map<const VoronoiDiagram::Vertex*, unsigned> vertex_reference{};
    std::list<VoronoiDiagram::Vertex> const& voronoi_vertices = voronoi_diagram.getVertices();

    for (auto const& vertex : voronoi_vertices)
    {
        vertex_indices[&vertex] = vertex_idx;
        ++vertex_idx;
    }

    std::vector<Node<2>*> nodes{};
    std::vector<VertexElement<2, 2>*> elements{};

    for (auto const& pair : vertex_indices)
    {
        auto const& position = pair.first->point;
        auto const node_idx = pair.second;
        Node<2>* p_node = new Node<2>{ node_idx, false, position.x, position.y, 0 };
        nodes.push_back(p_node);
    }

    for (unsigned i = 0; i < voronoi_diagram.getNbSites(); ++i)
    {
        std::vector<Node<2>*> element_nodes{};
        VoronoiDiagram::Face* face = voronoi_diagram.getFace(i);
        VoronoiDiagram::HalfEdge* half_edge = face->outerComponent;
        if (half_edge == nullptr)
        {
            continue;
        }
        while (half_edge->prev != nullptr)
        {
            half_edge = half_edge->prev;
            if (half_edge == face->outerComponent)
            {
                break;
            }
        }
        VoronoiDiagram::HalfEdge* const start = half_edge;
        while (half_edge != nullptr)
        {
            VoronoiDiagram::Vertex* origin = half_edge->origin;
            unsigned const node_idx = vertex_indices.at(origin);
            element_nodes.push_back(nodes[node_idx]);
            vertex_reference[origin] += 1;
            half_edge = half_edge->next;
            if (half_edge == start)
            {
                break;
            }
        }
        VertexElement<2, 2>* p_element = new VertexElement<2, 2>{ i, element_nodes };
        elements.push_back(p_element);
    }

    unsigned node_i = 0;
    for (auto const& pair : vertex_indices)
    {
        bool const on_boundary = vertex_reference.at(pair.first) <= 2;
        nodes[node_i]->SetAsBoundaryNode(on_boundary);
        ++node_i;
    }

    mpMesh = new MutableVertexMesh<2, 2>{ nodes, elements };
}

MyVoronoiVertexMeshGenerator::~MyVoronoiVertexMeshGenerator()
{
    delete mpMesh;
}

MutableVertexMesh<2, 2>* MyVoronoiVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}
