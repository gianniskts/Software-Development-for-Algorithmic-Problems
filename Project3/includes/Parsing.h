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
    int L = 50;
    double alpha = 2.0;
    double beta = 0.5;
    double kappa = 200.0;
    double xi = 1.1;
    double psi = 12.0;
    double lambda = 0.999;
    int random_deadlock_threshold = 8;  // # of no-improvement steps
    std::vector<std::string> steiner_methods;
    bool randomization_enabled = true;
};

// Parse input file in JSON format
InputJSON parse_file(const std::string& filename);

// Set the results to JSON format
void output_results(const std::string& filename, const InputJSON& input, const Triangulation& triangulation);

std::vector<std::string> split_csv(const std::string &s);

std::string vector_to_string(const std::vector<std::string>& vec);

#endif // PARSING_H
