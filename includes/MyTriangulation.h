#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include "../includes/PolygonManipulation.h"

// Class to hold and access triangulation info
class Triangulation {
public:
    CDT cdt;
    Polygon_2 polygon;
    std::unordered_map<CDT::Face_handle, bool> in_domain_map;
    boost::associative_property_map<std::unordered_map<CDT::Face_handle, bool>> in_domain;

    // Constructor
    Triangulation(CDT& cdt, Polygon_2& polygon) : cdt(cdt), polygon(polygon), in_domain(in_domain_map) {
    }

    // Copy constructor
    Triangulation(const Triangulation& other)
        : cdt(other.cdt),  // Copy the CDT
          polygon(other.polygon),  // Copy the Polygon
          in_domain_map(other.in_domain_map),  // Copy the in_domain_map
          in_domain(in_domain_map) {  // Reinitialize the property map
        
    }

    // Function to mark triangulation domain
    void mark_domain();

    // Function to check if a triangle is inside the boundary
    bool is_face_in_domain(CDT::Face_handle face) const;

    // Function to count triangulation's obtuse triangles
    int count_obtuse_triangles();
};

// Forward declaration
bool is_obtuse(const Point& A, const Point& B, const Point& C);

#endif // TRIANGULATION_H