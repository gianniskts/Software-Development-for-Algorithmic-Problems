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
        InputJSON<int> input_json = parse_file<int>(json_file);

        // Perform constrained delaunay triangulation
        result = delaunay_const_triangulation<int>(input_json);

        // Output the results to a JSON file
        output_results<int>("output.json", input_json, result);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
