#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include "../includes/PolygonManipulation.h"

// Class to hold and access triangulation info
class Triangulation : public CDT {
public:
    int min_obtuse_triangles;
    CDT cdt;
    Polygon_2 polygon;
    std::unordered_map<CDT::Face_handle, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<CDT::Face_handle, bool>> in_domain;
    bool randomizationUsed = false;
    double p_bar = 0.0;  

    // Default constructor
    Triangulation()
        : cdt(), polygon(), in_domain(in_domain_map) {
        randomizationUsed = false;
    }

    // Constructor with parameters
    Triangulation(CDT& cdt, Polygon_2& polygon)
        : cdt(cdt), polygon(polygon), in_domain(in_domain_map) {
        randomizationUsed = false;
    }

    // Copy constructor
    Triangulation(const Triangulation& other)
        : cdt(other.cdt),  // Copy the CDT
          polygon(other.polygon),  // Copy the Polygon
          in_domain_map(other.in_domain_map),  // Copy the in_domain_map
          in_domain(other.in_domain), // Copy the in_domain
          min_obtuse_triangles(other.min_obtuse_triangles), // Copy the min_obtuse_triangles
          randomizationUsed(other.randomizationUsed), // Copy the randomizationUsed
          p_bar(other.p_bar) { // Copy the p_bar
    }

    // Function to mark triangulation domain
    void mark_domain();

    // Function to check if a triangle is inside the boundary
    bool is_face_in_domain(CDT::Face_handle face) const;

    // Function to check if an edge is inside the boundary
    bool is_edge_in_domain(CDT::Edge edge) const;

    // Function to count triangulation's obtuse triangles
    int count_obtuse_triangles() const;
};


#endif // TRIANGULATION_H