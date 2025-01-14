#ifndef ACTION_FUNCTIONS_H
#define ACTION_FUNCTIONS_H

#include "../includes/MyTriangulation.h"
#include "../includes/PolygonManipulation.h"

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

// Function to try edge flips on cdt
void edge_flip(Triangulation& triangulation);

// Function to project point A onto the line defined by B and C
Point project_point_onto_line(const Point& A, const Point& B, const Point& C);

// Function to check if an edge is bound
bool is_edge_valid(const Point& v1, const Point& v2, const Polygon_2& polygon);

// Function to find the midpoint of the largest edge formed by the points given
Point get_midpoint(const Point& A, const Point& B, const Point& C);

// Function to add Steiner points
bool add_optimal_steiner(Triangulation& triangulation);

// Function to eliminate obtuse triangles (if possible)
void eliminate_obtuse_triangles(Triangulation& triangulation);

// Function to generate all subsets of a given set
std::vector<std::vector<std::string>> generate_subsets(const std::vector<std::string>& items);

// Function to detect the category of the input
std::string detect_category(const InputJSON& input);

// Function to insert a Steiner point near the centroid of an obtuse triangle with a Gaussian distribution (if possible)
bool insert_gaussian_near_centroid_of_obtuse_triangle(Triangulation &triang, double sigma);

// Function to compute p^n
double compute_pn(int obtuse_prev, int obtuse_curr, int n_prev, int n_curr);

// Function to randomize the triangulation if stuck in a local minimum for too long (if possible)
void randomize_if_stuck(Triangulation &triang, int stepsWithoutImprovement, int threshold, bool randomization_enabled);

#endif // ACTION_FUNCTIONS_H