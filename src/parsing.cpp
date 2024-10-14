#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/json/src.hpp>

#include "../includes/parsing.h"

#define CONTENT_TYPE "CG_SHOP_2025_Solution"

namespace json = boost::json;

// Function to parse JSON file
template <typename T>
InputJSON<T> parse_file(const std::string& filename) {
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
    InputJSON<T> data;
    data.instance_uid = json::value_to<std::string>(json_obj["instance_uid"]);
    data.num_points = json::value_to<int>(json_obj["num_points"]);

    // Parse the arrays
    for (const auto& val : json_obj["points_x"].as_array()) {
        data.points_x.push_back(json::value_to<T>(val));
    }

    for (const auto& val : json_obj["points_y"].as_array()) {
        data.points_y.push_back(json::value_to<T>(val));
    }

    for (const auto& val : json_obj["region_boundary"].as_array()) {
        data.region_boundary.push_back(json::value_to<T>(val));
    }

    data.num_constraints = json::value_to<int>(json_obj["num_constraints"]);

    // Parse additional_constraints, which is a list of lists
    for (const auto& constraint : json_obj["additional_constraints"].as_array()) {
        auto pair = constraint.as_array();
        data.additional_constraints.emplace_back(json::value_to<T>(pair[0]), json::value_to<T>(pair[1]));
    }

    return data;
}

// Function for results output
template <typename T>
void output_results(const std::string& filename, const InputJSON<T>& input, const Polygon_2& polygon) {
    // Object to store output json
    boost::json::object results;
    
    // Arrays to store results from the polygon
    json::array steiner_points_x;
    json::array steiner_points_y;
    json::array edges;

    // Set the values to corresponding fields
    results["content_type"] = CONTENT_TYPE;
    results["instance_uid"] = input.instance_uid;

    // Iterate any added point to the initial polygon
    for (size_t i = input.num_points; i < polygon.size(); ++i) {

        // Get exact coordinates
        const auto point_x = CGAL::exact(polygon[i].x());
        const auto point_y = CGAL::exact(polygon[i].y());

        std::ostringstream ss;

        // Check if the number is integer or rational and process accordingly
        if (point_x.get_den() == 1) {
            ss << point_x;
            steiner_points_x.emplace_back(stoi(ss.str()));
        } else {
            ss << point_x;
            steiner_points_x.emplace_back(ss.str());
        }
        ss.str("");


        if (point_y.get_den() == 1) {
            ss << point_y;
            steiner_points_y.emplace_back(stoi(ss.str()));
        } else {
            ss << point_y;
            steiner_points_y.emplace_back(ss.str());
        }
        ss.str("");
        
    }

    // Set the values to the corresponding fields
    results["steiner_points_x"] = steiner_points_x;
    results["steiner_points_y"] = steiner_points_y;

    
    /*
    for (const auto& constraint : input.additional_constraints) {
        json::array edge;
        edge.push_back(json::value(constraint.first));  // First point index
        edge.push_back(json::value(constraint.second)); // Second point index
        edges.push_back(edge);
    }*/

    // Set the values to the corresponding field
    results["edges"] = edges;

    // Output results
    std::string json_string = json::serialize(results);
    std::ofstream json_file(filename);

    // Security check if the file opens
    if (json_file.is_open()) {
        json_file << json_string;
        json_file.close();
    } else {
        throw std::runtime_error("Unable to open output JSON file for writing");
    }
    
}

template InputJSON<int> parse_file<int>(const std::string& filename);
template void output_results<int>(const std::string& filename, const InputJSON<int>& input, const Polygon_2& polygon);
