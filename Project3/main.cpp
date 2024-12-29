#include <iostream>
#include <string>
#include <vector>
#include "../includes/Parsing.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include "../includes/OptimizationMethods.h"

int main(int argc, char* argv[]) {
    try {
        if (argc != 5) {
            std::cerr << "Usage: ./build/main -i /path/to/input.json -o /path/to/output.json" << std::endl;
            return 1;
        }

        std::string input_file;
        std::string output_file;

        // Command-line argument handling
        for (int i = 1; i < argc; i += 2) {
            std::string arg = argv[i];
            if (arg == "-i") {
                input_file = argv[i + 1];
            } else if (arg == "-o") {
                output_file = argv[i + 1];
            } else {
                std::cerr << "Unknown argument: " << arg << std::endl;
                return 1;
            }
        }

        // Parse the input file
        InputJSON input_json = parse_file(input_file);

        std::function<Triangulation(const InputJSON&)> method_function;
        if (input_json.method == "local") {
            method_function = local_search;
        } else if (input_json.method == "sa") {
            method_function = simulated_annealing;
        } else if (input_json.method == "ant") {
            method_function = ant_colony_optimization;
        } else {
            throw std::runtime_error("Unknown method specified");
        }

        Triangulation result = method_function(input_json);

        // Output the results to the specified JSON file
        output_results(output_file, input_json, result);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
