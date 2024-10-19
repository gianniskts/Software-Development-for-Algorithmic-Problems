#ifndef PARSING_H
#define PARSING_H

#include <string>
#include <vector>
#include "../includes/PolygonManipulation.h"

// Struct to capture input data
template <typename T>
struct InputJSON {
    std::string instance_uid;
    int num_points;
    std::vector<T> points_x;
    std::vector<T> points_y;
    std::vector<T> region_boundary;
    int num_constraints;
    std::vector<std::pair<T, T>> additional_constraints;
};

// Parse input file in JSON format
template <typename T>
InputJSON<T> parse_file(const std::string& filename);

// Set the results to JSON format
template <typename T>
void output_results(const std::string& filename, const InputJSON<T>& input, const Polygon_2& polygon);

#endif // PARSING_H
