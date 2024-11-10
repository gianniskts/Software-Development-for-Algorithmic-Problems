// include/ActionFunctions.h

#ifndef ACTION_FUNCTIONS_H
#define ACTION_FUNCTIONS_H

#include "Triangulation.h"
#include "PolygonManipulation.h"

// Structure to store triangulation states
struct TriangulationState {
    std::vector<Point> steiner_points;
    int obtuse_triangle_count;

    // Custom comparator to prioritize fewer obtuse triangles
    bool operator<(const TriangulationState& other) const {
        return obtuse_triangle_count > other.obtuse_triangle_count;  // Min-heap based on obtuse_triangle_count
    }
};

// Function to check if a triangle is obtuse
bool is_obtuse(const Point& A, const Point& B, const Point& C);

// Function to try edge flips on the triangulation
void edge_flip(Triangulation& triangulation);

// Function to project point A onto the line defined by B and C
Point project_point_onto_line(const Point& A, const Point& B, const Point& C);

// Function to check if an edge is valid (inside the polygon)
bool is_edge_valid(const Point& v1, const Point& v2, const Polygon_2& polygon);

// Function to find the midpoint of the largest edge formed by the points given
Point get_midpoint(const Point& A, const Point& B, const Point& C);

// Function to add Steiner points
bool add_optimal_steiner(Triangulation& triangulation);

// Function to eliminate obtuse triangles (if possible)
void eliminate_obtuse_triangles(Triangulation& triangulation);

#endif // ACTION_FUNCTIONS_H
