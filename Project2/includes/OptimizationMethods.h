#ifndef OPTIMIZATION_METHODS_H
#define OPTIMIZATION_METHODS_H

#include "Parsing.h"
#include "MyTriangulation.h"

// Function declarations
Triangulation local_search(const InputJSON& input);
Triangulation simulated_annealing(const InputJSON& input);
Triangulation ant_colony_optimization(const InputJSON& input);

#endif // OPTIMIZATION_METHODS_H
