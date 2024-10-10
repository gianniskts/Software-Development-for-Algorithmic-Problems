#ifndef EDGEFLIP_H
#define EDGEFLIP_H

#include "../includes/dtriangulation.h"

// Function to check if a triangle is obtuse
bool is_obtuse(const Point& A, const Point& B, const Point& C);

// Function to iterate through CDT and count obtuse triangles
int count_obtuse_triangles(CDT& cdt);

// Function to eliminate obtuse triangles (if possible)
void eliminate_obtuse_triangles(CDT& cdt, Polygon_2& steiner_points);

#endif // EDGEFLIP_H