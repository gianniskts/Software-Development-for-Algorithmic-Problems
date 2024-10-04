#include "../includes/dtriangulation.h"

bool is_point_in_boundary(const Point& point, const Polygon_2 polygon) {
    return polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE;
}


// Function to set instance region boundary
void set_boundary(CDT& cdt, const std::vector<Point>& boundary) {
    // Represent the boundary as a polygon to restrict the triangulation
    Polygon_2 polygon(boundary.begin(), boundary.end());
    
    for (size_t i = 0; i < boundary.size(); i++) {
        // Add constraint between consecutive boundary points
        cdt.insert_constraint(boundary[i], boundary[(i + 1) % boundary.size()]);
    }

    // TODO: Set region boundary as a polygon
    // in order to restrict the triangulation inside it

    // Extra TODO: na spasw ton kwdika se perissotera files

}

// Function to perform delaunay triangulation
void delaunay_const_triangulation(InputJSON input_data) {
    
    // Create CDT
    CDT cdt;

    // Store input data
    std::vector<Point> points;
    std::vector<Point> region_boundary;

    for (size_t i = 0; i < input_data.num_points; i++) {
        points.push_back(Point(input_data.points_x[i], input_data.points_y[i]));
    }

    for (size_t i = 0; i < input_data.region_boundary.size(); i++) {
        region_boundary.push_back(points[input_data.region_boundary[i]]);
    }

    set_boundary(cdt, region_boundary);

    for (size_t i = 0; i < input_data.num_points; i++) {
        cdt.insert(points[i]);
    }

    CGAL::make_conforming_Delaunay_2(cdt);

    // Visualize CDT's results
    CGAL::draw(cdt);
}