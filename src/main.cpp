#include <iostream>
#include <string>
#include <vector>
#include "../includes/parsing.h"
#include "../includes/dtriangulation.h"

int main() {
    try {
        // Parse the input file
        InputJSON input_json = parse_file("../cgshop2025_examples_ortho_10_ff68423e.json");

        // Perform constrained delaunay triangulation
        delaunay_const_triangulation(input_json);

        // So far only prints back input json (to be fixed after steiner)
        output_results("../output.json", input_json);

    } 
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}