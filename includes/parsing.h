#ifndef PARSING_H
#define PARSING_H

#include <string>
#include <vector>

// Struct to capture input data
struct InputJSON {
    std::string instance_uid;
    int num_points;
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<int> region_boundary;
    int num_constraints;
    std::vector<std::pair<int, int>> additional_constraints;
};

// Parse input file in JSON format
InputJSON parse_file (const std::string& filename);

// Output project's results in JSON format
void output_results (const std::string& filename, const InputJSON& input);

#endif // PARSING_H