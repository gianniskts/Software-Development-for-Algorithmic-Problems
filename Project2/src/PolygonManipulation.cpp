#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"
#include "../includes/Parsing.h"
#include "draw.h"


// Function to perform delaunay triangulation and visualize results
Triangulation delaunay_const_triangulation(const InputJSON input_data) {
    
    CDT cdt;
    Polygon_2 polygon;
    std::vector<Point> region_boundary;
    std::vector<std::pair<Point, Point>> extra_constraints;
    
    // Store polygon's points
    for (size_t i = 0; i < input_data.num_points; ++i) {
        polygon.push_back(Point(input_data.points_x[i], input_data.points_y[i]));
        cdt.insert(polygon[i]);
    }

    // Store boundary points
    for (size_t i = 0; i < input_data.region_boundary.size(); ++i) {
        region_boundary.push_back(polygon[input_data.region_boundary[i]]);
    }

    // Set region boundary
    cdt.insert_constraint(region_boundary.begin(), region_boundary.end(), true);

    // If additional contraints are given
    if (input_data.num_constraints) {
        for (size_t i = 0; i < input_data.num_constraints; ++i) {
            Point p1 = polygon[input_data.additional_constraints[i].first];
            Point p2 = polygon[input_data.additional_constraints[i].second];
            extra_constraints.push_back(std::make_pair(p1, p2));
        }

        // Set extra constraints
        for (size_t i = 0; i < extra_constraints.size(); ++i) {
            cdt.insert_constraint(extra_constraints[i].first, extra_constraints[i].second);
        }
    }

    // Construct a triangulation instance with the given parameters
    Triangulation triangulation(cdt, polygon);
    
    //Mark facets that are inside the domain bounded by the polygon
    triangulation.mark_domain();

    triangulation.min_obtuse_triangles = triangulation.count_obtuse_triangles();
    
    // Try edge flips on obtuse triangles
    edge_flip(triangulation);

    // Attempt to eliminate obtuse triangles by adding Steiner points
    eliminate_obtuse_triangles(triangulation, input_data.L);
    
    // Visualize CDT's results
    CGAL::draw(triangulation.cdt, triangulation.in_domain);
    

    return triangulation;
}