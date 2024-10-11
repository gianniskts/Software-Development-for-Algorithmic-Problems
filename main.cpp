#include <iostream>
#include <string>
#include <vector>
#include "../includes/parsing.h"
#include "../includes/dtriangulation.h"
#include "../includes/edgeflip.h"

int main() {
    try {
        // Parse the input file
        InputJSON<int> input_json = parse_file<int>("../cgshop2025_examples_ortho_10_ff68423e.json");

        // Perform constrained delaunay triangulation
        delaunay_const_triangulation<int>(input_json);

        // Output the results to a JSON file
        output_results<int>("../output.json", input_json);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
