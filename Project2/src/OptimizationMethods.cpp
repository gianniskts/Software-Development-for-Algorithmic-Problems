#include "../includes/OptimizationMethods.h"
#include "../includes/PolygonManipulation.h"

Triangulation local_search(const InputJSON& input) {
    return delaunay_const_triangulation(input);
}

Triangulation simulated_annealing(const InputJSON& input) {
    // For now, returning a basic triangulation
    return delaunay_const_triangulation(input);
}

Triangulation ant_colony_optimization(const InputJSON& input) {
    // For now, returning a basic triangulation
    return delaunay_const_triangulation(input);
}
