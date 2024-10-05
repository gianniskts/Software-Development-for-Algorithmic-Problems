// TODO: Search how to add Steiner points !!
// TODO: Fix output json after Steiner pts ~

#include "../includes/dtriangulation.h"
#include "../includes/rboundary.h"
#include "../includes/xtrconstraints.h"

// Function to perform delaunay triangulation and visualize results
void delaunay_const_triangulation(InputJSON input_data) {
    
    // Create CDT
    CDT cdt;

    // Store input data
    std::vector<Point> points;
    std::vector<Point> region_boundary;
    std::vector<std::pair<Point, Point>> extra_constraints;

    for (size_t i = 0; i < input_data.num_points; i++) {
        points.push_back(Point(input_data.points_x[i], input_data.points_y[i]));
    }

    // If region boundary is given
    if (input_data.region_boundary.size() != 0) {
        for (size_t i = 0; i < input_data.region_boundary.size(); i++) {
            region_boundary.push_back(points[input_data.region_boundary[i]]);
        }

        // Set region boundary
        set_boundary(cdt, region_boundary);
    }

    // If additional contraints are given
    if (input_data.num_constraints) {
        for (size_t i = 0; i < input_data.num_constraints; i++) {
            Point p1 = points[input_data.additional_constraints[i].first];
            Point p2 = points[input_data.additional_constraints[i].second];
            extra_constraints.push_back(std::make_pair(p1, p2));
        }

        // Set additional constraints
        set_constraints(cdt, extra_constraints);
    }    

    // Add all the points in the triangulation
    for (size_t i = 0; i < input_data.num_points; i++) {
        cdt.insert(points[i]);
    }

    // Refines the constrained DT into a conforming DT
    CGAL::make_conforming_Delaunay_2(cdt);

    // Visualize CDT's results
    CGAL::draw(cdt);
}