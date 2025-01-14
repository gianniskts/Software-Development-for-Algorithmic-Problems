#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/json/src.hpp>
#include "../includes/ActionFunctions.h"
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

        if (data.parameters.count("L") > 0) {
            data.L = json::value_to<int>(data.parameters["L"]);
        }

        if (data.parameters.count("alpha") > 0) {
            data.alpha = json::value_to<double>(data.parameters["alpha"]);
        }

        if (data.parameters.count("beta") > 0) {
            data.beta = json::value_to<double>(data.parameters["beta"]);
        }
        if (data.parameters.count("xi") > 0) {
            data.xi = json::value_to<double>(data.parameters["xi"]);
        }
        if (data.parameters.count("psi") > 0) {
            data.psi = json::value_to<double>(data.parameters["psi"]);
        }
        if (data.parameters.count("lambda") > 0) {
            data.lambda = json::value_to<double>(data.parameters["lambda"]);
        }
        if (data.parameters.count("kappa") > 0) {
            data.kappa = json::value_to<int>(data.parameters["kappa"]);
        }
        // Possibly parse a random_deadlock_threshold
        if (data.parameters.count("random_deadlock_threshold") > 0) {
            data.random_deadlock_threshold =
                json::value_to<int>(data.parameters["random_deadlock_threshold"]);
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
void output_results(const std::string& filename, const InputJSON& input, const Triangulation& triangulation, bool advancedOutput) {
    
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
    params_json["alpha"]  = input.alpha;
    params_json["beta"]   = input.beta;
    params_json["L"]      = input.L;
    params_json["xi"]     = input.xi;
    params_json["psi"]    = input.psi;
    params_json["lambda"] = input.lambda;
    params_json["kappa"]  = input.kappa;

    // Overwrite with any user-supplied parameters from input.parameters
    // (This merges any JSON parse results with the final numeric values.)
    for (const auto &p : input.parameters) {
        params_json[p.first] = p.second;
    }

    // Now store it
    results["parameters"] = params_json;
    results["randomization"] = triangulation.randomizationUsed;

    int num_steiner_points = triangulation.cdt.number_of_vertices() - input.num_points;
    double final_energy = input.alpha * obtuse_count + input.beta * num_steiner_points;
    std::ostringstream energy_stream;
    energy_stream << std::fixed << std::setprecision(2) << final_energy;
    if(advancedOutput)
        results["energy"] = energy_stream.str();
    std::cout << "[DEBUG] Final Energy: " << final_energy << " | Obtuse Count: " 
          << obtuse_count << " | Steiner Points: " << num_steiner_points << std::endl;
    if(advancedOutput)
        results["p_bar"] = triangulation.p_bar;

    std::string category = detect_category(input);
    if(advancedOutput)
        results["category"] = category;

    // Output results
    file << json::serialize(results);
    file.close();
}

std::vector<std::string> split_csv(const std::string &s) {
    // Helper: splits a comma-separated string (e.g. "local,sa,ant") into tokens
    std::vector<std::string> result;
    std::stringstream ss(s);
    while (ss.good()) {
        std::string substr;
        if (!std::getline(ss, substr, ',')) break;
        // trim spaces
        substr.erase(0, substr.find_first_not_of(" \t\r\n"));
        substr.erase(substr.find_last_not_of(" \t\r\n") + 1);
        if (!substr.empty()) result.push_back(substr);
    }
    return result;
}

std::string vector_to_string(const std::vector<std::string>& vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i < vec.size() - 1) {
            oss << ", ";
        }
    }
    oss << "]";
    return oss.str();
}

