#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/json/src.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "../includes/parsing.h"

#define CONTENT_TYPE "CG_SHOP_2025_Solution"

namespace json = boost::json;

// Function to parse JSON file
InputJSON parse_file (const std::string& filename) {
    
    // Open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open JSON file");
    }

    // Parse the JSON file as a string
    std::ostringstream ss;
    ss << file.rdbuf();
    std::string json_content = ss.str();

    json::value json_value = json::parse(json_content);
    json::object json_obj = json_value.as_object();

    // Populate the struct
    InputJSON data;
    data.instance_uid = json::value_to<std::string>(json_obj["instance_uid"]);
    data.num_points = json::value_to<int>(json_obj["num_points"]);

    // Parse the arrays
    for (const auto& val : json_obj["points_x"].as_array()) {
        data.points_x.push_back(json::value_to<int>(val));
    }

    for (const auto& val : json_obj["points_y"].as_array()) {
        data.points_y.push_back(json::value_to<int>(val));
    }

    for (const auto& val : json_obj["region_boundary"].as_array()) {
        data.region_boundary.push_back(json::value_to<int>(val));
    }

    data.num_constraints = json::value_to<int>(json_obj["num_constraints"]);

    // Parse additional_constraints, which is a list of lists
    for (const auto& constraint : json_obj["additional_constraints"].as_array()) {
        auto pair = constraint.as_array();
        data.additional_constraints.emplace_back(json::value_to<int>(pair[0]), json::value_to<int>(pair[1]));
    }

    return data;
}

//Function for results output
void output_results (const std::string& filename, const InputJSON& input) {
    boost::property_tree::ptree results;

    results.put("content_type", CONTENT_TYPE);
    results.put("instance_uid", input.instance_uid);

    // example
    json::array steiner_points_x;
    json::array steiner_points_y;

    // For demonstration, I just copy the input points as steiner points
    for (const auto& x : input.points_x) {
        steiner_points_x.push_back(json::value(x));  // Integer points
    }
    for (const auto& y : input.points_y) {
        steiner_points_y.push_back(json::value(y));  // Integer points
    }

    results.put("steiner_points_x", steiner_points_x);
    results.put("steiner_points_y", steiner_points_y);

    json::array edges;
    for (const auto& constraint : input.additional_constraints) {
        json::array edge;
        edge.push_back(json::value(constraint.first));  // First point index
        edge.push_back(json::value(constraint.second)); // Second point index
        edges.push_back(edge);
    }
    results.put("edges", edges);

    // Set results on output json file
    std::ofstream output_file(filename);
    if (!output_file.is_open()) {
        throw std::runtime_error("Unable to open output JSON file for writing");
    }
    boost::property_tree::write_json(output_file, results);
    output_file.close();
}
