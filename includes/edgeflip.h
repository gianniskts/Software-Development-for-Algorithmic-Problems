#ifndef EDGEFLIP_H
#define EDGEFLIP_H

#include "../includes/dtriangulation.h"

// Function to calculate angle at a vertex
double angle(const Point& p1, const Point& p2, const Point& p3);

// Function to check if we must flip edge
bool should_flip(const CDT& cdt, Edge edge);

bool is_flip_valid(const CDT& cdt, const CDT::Edge& edge);

// Function to perform edge flips
void edge_flip(CDT& cdt);

#endif // EDGEFLIP_H