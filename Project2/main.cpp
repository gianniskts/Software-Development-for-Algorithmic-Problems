// src/main.cpp

#include <iostream>
#include <string>
#include <vector>
#include "../include/Parsing.h"
#include "../include/PolygonManipulation.h"
#include "../include/ActionFunctions.h"
#include "../include/Triangulation.h"

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cerr << "Usage: ./program <input_json_file>" << std::endl;
            return 1;
        }

        std::string json_file = argv[1];

        // Parse the input file
        InputJSON input_json = parse_file(json_file);

        // Perform constrained delaunay triangulation and optimization
        Triangulation result = delaunay_const_triangulation(input_json);

        // Output the results to a JSON file
        output_results("output.json", input_json, result);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
