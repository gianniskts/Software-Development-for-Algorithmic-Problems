#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/json/src.hpp>

#include "../includes/Parsing.h"

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

    // Parse the method
    if (json_obj.contains("method")) {
        data.method = json::value_to<std::string>(json_obj["method"]);
    } else {
        data.method = "local";
    }

    // Parse the parameters
    if (json_obj.contains("parameters")) {
        for (const auto& param : json_obj["parameters"].as_object()) {
            data.parameters[param.key_c_str()] = param.value();
        }
    }

    // Parse the delaunay flag
    if (json_obj.contains("delaunay")) {
        data.delaunay = json::value_to<bool>(json_obj["delaunay"]);
    } else {
        data.delaunay = true;
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
    boost::json::object results;
    
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
        const auto point_x = CGAL::exact(points[i].x());
        const auto point_y = CGAL::exact(points[i].y());

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

    // Iterate through the edges of the result triangulation store the points indices
    for (auto eit = triangulation.cdt.edges_begin(); eit != triangulation.cdt.edges_end(); ++eit) {
        const Edge& edge = *eit;

        // Get all edges of triangulation (constraints as well)
        if (triangulation.is_edge_in_domain(edge) || triangulation.cdt.is_constrained(edge)) {
            auto source_it = std::find(points.begin(), points.end() - 1, edge.first->vertex((edge.second + 1) % 3)->point());
            auto target_it = std::find(points.begin(), points.end() - 1, edge.first->vertex((edge.second + 2) % 3)->point());

            int source_index = std::distance(points.begin(), source_it);
            int target_index = std::distance(points.begin(), target_it);

            edges.emplace_back(json::array{source_index, target_index});
        }
        
    }

    // Set the values to the corresponding field
    results["edges"] = edges;

    // Include the obtuse_count
    int obtuse_count = triangulation.count_obtuse_triangles();
    results["obtuse_count"] = obtuse_count;

    // Include method and parameters
    results["method"] = input.method;

    // Serialize parameters back to JSON
    json::object params_json;
    // TODO: L is an int and needs better handling
    for (const auto& param : input.parameters) {
        params_json[param.first] = param.second;
    }
    results["parameters"] = params_json;

    // Output results
    file << json::serialize(results);
    file.close();
}
