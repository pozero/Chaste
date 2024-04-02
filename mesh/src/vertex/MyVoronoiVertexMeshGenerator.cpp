#include "MyVoronoiVertexMeshGenerator.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_adaptation_policies_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Weighted_point_2.h>
#include <CGAL/bounding_box.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using RT = CGAL::Regular_triangulation_2<K>;
using AT = CGAL::Regular_triangulation_adaptation_traits_2<RT>;
using DRP = CGAL::Regular_triangulation_degeneracy_removal_policy_2<RT>;
using VD = CGAL::Voronoi_diagram_2<RT, AT, DRP>;
using vector = K::Vector_2;
using point = K::Point_2;
using weight_point = K::Weighted_point_2;

MyVoronoiVertexMeshGenerator::MyVoronoiVertexMeshGenerator(
    std::vector<double>& cell_target_areas)
{
    uint32_t const cell_count = cell_target_areas.size();
    RandomNumberGenerator* rng = RandomNumberGenerator::Instance();
    std::vector<point> sites{};
    std::vector<weight_point> weighted_sites{};
    std::map<point, int32_t> site_indices{};
    sites.reserve(cell_count);
    weighted_sites.reserve(cell_count);
    double total_area = 0.0;
    for (double const area : cell_target_areas)
    {
        total_area += area;
    }
    total_area *= 0.5;
    double const average_area = total_area / (double)cell_count;
    double const site_coord_range = std::sqrt(total_area);
    double const site_side_length = std::sqrt(average_area);
    double const site_coord_min = site_side_length;
    for (uint32_t si = 0; si < cell_count; ++si)
    {
        double const x = site_coord_range * rng->ranf() + site_coord_min;
        double const y = site_coord_range * rng->ranf() + site_coord_min;
        sites.push_back(point{ x, y });
        weighted_sites.push_back(weight_point{ point{ x, y }, cell_target_areas[si] });
        site_indices.emplace(point{ x, y }, (int32_t)si);
    }
    std::vector<point> vertices{};
    std::vector<bool> vertex_at_boundary{};
    std::map<point, uint32_t> vertex_indices{};
    std::vector<std::vector<uint32_t> > site_vertices(cell_count, std::vector<uint32_t>{});
    RT triangulation{ weighted_sites.begin(), weighted_sites.end() };
    if (!triangulation.is_valid())
    {
        EXCEPTION("Invalid triangulation\n");
    }
    auto const bounding_box = CGAL::bounding_box(sites.begin(), sites.end());
    double const box_side_length = std::max(bounding_box.xmax() - bounding_box.xmax(),
                                            bounding_box.ymax() - bounding_box.ymin());
    point center = CGAL::ORIGIN;
    for (uint32_t i = 0; i < weighted_sites.size(); ++i)
    {
        center += (sites[i] - CGAL::ORIGIN);
    }
    center = CGAL::ORIGIN + (center - CGAL::ORIGIN) / (double)sites.size();
    uint32_t const exterior_point_side_count = 8;
    double const exterior_point_radius = 0.65;
    std::vector<weight_point> exterior_points{};
    for (uint32_t i = 1; i <= exterior_point_side_count; ++i)
    {
        double const varying_coord = 2.0 * ((double)i / (double)exterior_point_side_count) * exterior_point_radius - exterior_point_radius;
        vector const displacement0{ exterior_point_radius, varying_coord };
        vector const displacement1{ -exterior_point_radius, varying_coord };
        vector const displacement2{ varying_coord, exterior_point_radius };
        vector const displacement3{ varying_coord, -exterior_point_radius };
        weight_point const exterior_p0{ center + box_side_length * displacement0,
                                        average_area };
        weight_point const exterior_p1{ center + box_side_length * displacement1,
                                        average_area };
        weight_point const exterior_p2{ center + box_side_length * displacement2,
                                        average_area };
        weight_point const exterior_p3{ center + box_side_length * displacement3,
                                        average_area };
        site_indices.insert(std::make_pair(exterior_p0.point(), -1));
        site_indices.insert(std::make_pair(exterior_p1.point(), -1));
        site_indices.insert(std::make_pair(exterior_p2.point(), -1));
        site_indices.insert(std::make_pair(exterior_p3.point(), -1));
        exterior_points.push_back(exterior_p0);
        exterior_points.push_back(exterior_p1);
        exterior_points.push_back(exterior_p2);
        exterior_points.push_back(exterior_p3);
    }
    triangulation.insert(exterior_points.begin(), exterior_points.end());
    if (!triangulation.is_valid())
    {
        EXCEPTION("Invalid triangulation\n");
    }
    for (auto face_iter = triangulation.faces_begin();
         face_iter != triangulation.faces_end(); ++face_iter)
    {
        point const& v0 = face_iter->vertex(0)->point().point();
        point const& v1 = face_iter->vertex(1)->point().point();
        point const& v2 = face_iter->vertex(2)->point().point();
        int32_t const s0 = site_indices.at(v0);
        int32_t const s1 = site_indices.at(v1);
        int32_t const s2 = site_indices.at(v2);
        if (s0 == -1 && s1 == -1 && s2 == -1)
        {
            continue;
        }
        bool const at_boundary = s0 == -1 || s1 == -1 || s2 == -1;
        point const dual_point = triangulation.dual(face_iter);
        uint32_t const new_vertex_idx = vertices.size();
        vertices.push_back(dual_point);
        vertex_at_boundary.push_back(at_boundary);
        vertex_indices.insert(std::make_pair(dual_point, new_vertex_idx));
    }
    VD voronoi{ triangulation };
    for (auto face_iter = voronoi.faces_begin(); face_iter != voronoi.faces_end();
         ++face_iter)
    {
        point const& site = face_iter->dual()->point().point();
        int32_t const site_idx = site_indices.at(site);
        if (site_idx < 0)
        {
            continue;
        }
        auto const face_half_begin = face_iter->halfedge();
        auto half = face_half_begin;
        do
        {
            point const& source_vertex = half->source()->point();
            uint32_t const vertex_idx = vertex_indices.at(source_vertex);
            site_vertices.at(site_idx).push_back(vertex_idx);
            half = half->next();
        } while (half != face_half_begin);
    }

    std::vector<Node<2>*> nodes{};
    std::vector<VertexElement<2, 2>*> elements{};
    for (uint32_t i = 0; i < vertices.size(); ++i)
    {
        point const& position = vertices[i];
        Node<2>* node = new Node<2>{ i, vertex_at_boundary[i], position.x(), position.y(), 0.0 };
        nodes.push_back(node);
    }
    std::vector<double> valid_cell_target_area{};
    for (uint32_t i = 0; i < sites.size(); ++i)
    {
        std::vector<uint32_t> const& site_vertex = site_vertices[i];
        if (site_vertex.empty())
        {
            continue;
        }
        valid_cell_target_area.push_back(weighted_sites[i].weight());
        std::vector<Node<2>*> element_nodes{};
        element_nodes.reserve(site_vertex.size());
        for (uint32_t const vertex_idx : site_vertex)
        {
            element_nodes.push_back(nodes[vertex_idx]);
        }
        VertexElement<2, 2>* element = new VertexElement<2, 2>{ (unsigned)elements.size(), element_nodes };
        elements.push_back(element);
    }
    mpMesh = new MutableVertexMesh<2, 2>{ nodes, elements };
    cell_target_areas = std::move(valid_cell_target_area);
}

MyVoronoiVertexMeshGenerator::~MyVoronoiVertexMeshGenerator()
{
    delete mpMesh;
}

MutableVertexMesh<2, 2>* MyVoronoiVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}
