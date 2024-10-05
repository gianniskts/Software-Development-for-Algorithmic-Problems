#include "../includes/rboundary.h"

// Function to set instance region boundary
void set_boundary(CDT& cdt, const std::vector<Point>& boundary) {
    // Represent the boundary as a polygon to restrict the triangulation
    Polygon_2 polygon(boundary.begin(), boundary.end());
    
    cdt.insert_constraint(polygon.begin(), polygon.end(), true);
}