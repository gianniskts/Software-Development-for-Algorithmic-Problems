#include "../includes/xtrconstraints.h"

// Function to set the additional constraints
void set_constraints(CDT& cdt, const std::vector<std::pair<Point,Point>>& constraints) {
    for (size_t i = 0; i < constraints.size(); i++) {
        cdt.insert_constraint(constraints[i].first, constraints[i].second);
    }
}