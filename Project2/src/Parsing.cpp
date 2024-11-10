// src/Parsing.cpp

#include "../include/Parsing.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/json/src.hpp>

#define CONTENT_TYPE "CG_SHOP_2025_Solution"

namespace json = boost::json;

// Function to parse JSON file
InputJSON parse_file(const std::string& filename) {
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

// Function for results output
void output_results(const std::string& filename, const InputJSON& input, const Triangulation& triangulation) {
    std::ofstream file(filename);

    // Check if file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "' for writing." << std::endl;
        return;
    }

    // Object to store output json
    json::object results;

    // Arrays to store results from the polygon
    json::array steiner_points_x;
    json::array steiner_points_y;
    json::array edges;

    // Set the values to corresponding fields
    results["content_type"] = CONTENT_TYPE;
    results["instance_uid"] = input.instance_uid;

    std::vector<Point> points;
    for (auto vit = triangulation.cdt.vertices_begin(); vit != triangulation.cdt.vertices_end(); ++vit) {
        points.push_back(vit->point());
    }

    // Iterate any added point to the initial polygon
    for (size_t i = input.num_points; i < points.size(); ++i) {
        // Get exact coordinates
        const auto point_x = CGAL::to_double(points[i].x());
        const auto point_y = CGAL::to_double(points[i].y());

        steiner_points_x.emplace_back(point_x);
        steiner_points_y.emplace_back(point_y);
    }

    // Set the values to the corresponding fields
    results["steiner_points_x"] = steiner_points_x;
    results["steiner_points_y"] = steiner_points_y;

    // Map points to indices
    std::map<Point, int> point_indices;
    int index = 0;
    for (const auto& pt : points) {
        point_indices[pt] = index++;
    }

    // Iterate through the edges of the result triangulation store the points indices
    for (auto eit = triangulation.cdt.edges_begin(); eit != triangulation.cdt.edges_end(); ++eit) {
        const Edge& edge = *eit;

        // Get all edges of triangulation (constraints as well)
        if (triangulation.is_edge_in_domain(edge) || triangulation.cdt.is_constrained(edge)) {
            Point source = edge.first->vertex((edge.second + 1) % 3)->point();
            Point target = edge.first->vertex((edge.second + 2) % 3)->point();

            int source_index = point_indices[source];
            int target_index = point_indices[target];

            edges.emplace_back(json::array{source_index, target_index});
        }
    }

    // Set the values to the corresponding field
    results["edges"] = edges;

    // Output results
    std::string json_string = json::serialize(results);

    // Security check if the file opens
    if (file.is_open()) {
        file << json_string;
        file.close();
    } else {
        throw std::runtime_error("Unable to open output JSON file for writing");
    }
}
