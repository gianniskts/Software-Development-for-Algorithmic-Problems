#ifndef PARSING_H
#define PARSING_H

#include <string>
#include <vector>
#include <map>
#include <boost/json.hpp>
#include "../includes/PolygonManipulation.h"
#include "../includes/MyTriangulation.h"

// Struct to capture input data
struct InputJSON {
    std::string instance_uid;
    int num_points;
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<int> region_boundary;
    int num_constraints;
    std::vector<std::pair<int, int>> additional_constraints;
    std::string method;
    std::map<std::string, boost::json::value> parameters;
    bool delaunay;
    int L;
    double alpha = 1.0;
    double beta = 1.0;
    double kappa = 10.0;
    double xi = 1.0;
    double psi = 1.0;
    double lambda = 1.0;
};

// Parse input file in JSON format
InputJSON parse_file(const std::string& filename);

// Set the results to JSON format
void output_results(const std::string& filename, const InputJSON& input, const Triangulation& triangulation);

#endif // PARSING_H
