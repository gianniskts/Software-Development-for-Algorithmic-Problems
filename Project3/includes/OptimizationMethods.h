#ifndef OPTIMIZATION_METHODS_H
#define OPTIMIZATION_METHODS_H

#include "Parsing.h"
#include "MyTriangulation.h"

struct AutoMethodResult {
    Triangulation triang;
    InputJSON best_input; 
};

// Function declarations
Triangulation local_search(const InputJSON& input);
Triangulation simulated_annealing(const InputJSON& input);
Triangulation ant_colony_optimization(const InputJSON& input);
AutoMethodResult auto_method(const InputJSON& input);

#endif // OPTIMIZATION_METHODS_H
