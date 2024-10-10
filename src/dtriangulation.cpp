// TODO: Mark anything outside the polygon as infinite
// Rel: Check example 8.4 (you might need to redefine a class for the graphics)
// TODO: Search how to add Steiner points !!
// TODO: Fix output json after Steiner pts ~
// TODO: Implement the algorithm using Templates !! (important)

#include "../includes/dtriangulation.h"
#include "../includes/edgeflip.h"
#include "draw.h"


// Function to perform delaunay triangulation and visualize results
void delaunay_const_triangulation(InputJSON input_data) {
    
    // Create CDT
    CDT cdt;

    // Store input data
    Polygon_2 polygon;
    std::vector<Point> region_boundary;
    std::vector<std::pair<Point, Point>> extra_constraints;
    

    // Store polygon's points
    for (size_t i = 0; i < input_data.num_points; i++) {
        polygon.push_back(Point(input_data.points_x[i], input_data.points_y[i]));

        cdt.insert(polygon[i]);
        
    }

    // Store boundary points
    for (size_t i = 0; i < input_data.region_boundary.size(); i++) {
        region_boundary.push_back(polygon[input_data.region_boundary[i]]);
    }

    // Set region boundary
    cdt.insert_constraint(region_boundary.begin(), region_boundary.end(), true);

    std::unordered_map<Face_handle, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<Face_handle,bool>> in_domain(in_domain_map);

    //Mark facets that are inside the domain bounded by the polygon
    


    // If additional contraints are given
    if (input_data.num_constraints) {
        for (size_t i = 0; i < input_data.num_constraints; i++) {
            Point p1 = polygon[input_data.additional_constraints[i].first];
            Point p2 = polygon[input_data.additional_constraints[i].second];
            extra_constraints.push_back(std::make_pair(p1, p2));
        }

        // Set extra constraints
        for (size_t i = 0; i < extra_constraints.size(); i++) {
            cdt.insert_constraint(extra_constraints[i].first, extra_constraints[i].second);
        }
    }

    eliminate_obtuse_triangles(cdt, polygon);

    CGAL::mark_domain_in_triangulation(cdt, in_domain);

    // Visualize CDT's results
    CGAL::draw(cdt, in_domain);
}