#include <iostream>
#include <string>
#include <vector>
#include "../include/parsing.h"

int main() {
    try {
        // Parse the input file
        InputJSON input_json = parse_file("input.json");

        output_results("output.json", input_json);

    } 
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}