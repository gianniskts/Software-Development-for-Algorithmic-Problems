#include <iostream>
#include <string>
#include <vector>
#include "../includes/Parsing.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"

int main(int argc, char* argv[]) {
    try {

        Polygon_2 result;

        std::string json_file = argv[1];

        // Parse the input file
        InputJSON input_json = parse_file(json_file);

        // Perform constrained delaunay triangulation
        result = delaunay_const_triangulation(input_json);

        // Output the results to a JSON file
        output_results("output.json", input_json, result);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
