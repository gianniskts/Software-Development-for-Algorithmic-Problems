#include <iostream>
#include <string>
#include <vector>
#include "../includes/parsing.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"

int main() {
    try {

        Polygon_2 result;

        // Parse the input file
        InputJSON<int> input_json = parse_file<int>("../input.json");

        // Perform constrained delaunay triangulation
        result = delaunay_const_triangulation<int>(input_json);

        // Output the results to a JSON file
        output_results<int>("../output.json", input_json, result);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
