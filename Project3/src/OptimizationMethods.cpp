#include "../includes/OptimizationMethods.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include <random>
#include <chrono>
#include <cmath>

#include "../includes/draw.h"

Triangulation local_search(const InputJSON& input) {
    Triangulation triangulation = delaunay_const_triangulation(input);
    
    bool improved = true;
    int iterations = 0;
    int L = input.L;

    while (improved && (iterations < L)) {
        //Mark facets that are inside the domain bounded by the polygon
        triangulation.mark_domain();

        // Flag to indicate if the triangulation was improved
        bool improved = false;

        Point best_steiner_point;
        int min_obtuse_triangles = triangulation.min_obtuse_triangles;   

        // Store the initial obtuse triangles
        std::vector<Face_handle> obtuse_faces;
        for (auto face_it = triangulation.cdt.finite_faces_begin(); face_it != triangulation.cdt.finite_faces_end(); ++face_it) {
            if (triangulation.is_face_in_domain(face_it) && is_obtuse(face_it->vertex(0)->point(), face_it->vertex(1)->point(), face_it->vertex(2)->point())) {
                obtuse_faces.push_back(face_it);
            }
        }

        // Iterate through the obtuse triangles
        for (const auto& face : obtuse_faces) {

            // Retrieve the triangle's points
            Face_handle fh = face;
            Point p0 = fh->vertex(0)->point();
            Point p1 = fh->vertex(1)->point();
            Point p2 = fh->vertex(2)->point();
            
            Point obtuse_angle, edge_vertex_1, edge_vertex_2;

            // Mark which angle is the obtuse
            if (angle(p0, p1, p2) == CGAL::OBTUSE) {
                obtuse_angle = p1;
                edge_vertex_1 = p2;
                edge_vertex_2 = p0;
            } else if (angle(p2, p0, p1) == CGAL::OBTUSE) {
                obtuse_angle = p0;
                edge_vertex_1 = p1;
                edge_vertex_2 = p2;
            } else {
                obtuse_angle = p2;
                edge_vertex_1 = p0;
                edge_vertex_2 = p1;
            }

            // Get multiple candidate steiner points positions
            Triangle_2 triangle(obtuse_angle, edge_vertex_1, edge_vertex_2);
            Point circumcenter = CGAL::circumcenter(triangle);
            Point centroid = CGAL::centroid(triangle);
            Point projection = project_point_onto_line(obtuse_angle, edge_vertex_1, edge_vertex_2);
            Point midpoint = get_midpoint(p0, p1, p2);
            Point polygon_centroid;

            // Vector to store candidate Steiner points
            std::vector<Point> candidate_points;

            // Validation checks to stay inbound
            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                && !triangulation.cdt.is_infinite(triangulation.cdt.locate(midpoint))) {
                candidate_points.push_back(midpoint);
            }
            
            if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(circumcenter))) {
                candidate_points.push_back(circumcenter);
            }
            
            if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(centroid))) {
                candidate_points.push_back(centroid);
            }

            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                && !triangulation.cdt.is_infinite(triangulation.cdt.locate(projection))) {
                candidate_points.push_back(projection);
            }

            // Iterate through the candidate points
            for (Point& candidate : candidate_points) {
                // Create a copy of the triangulation for testing
                Triangulation tcopy(triangulation);
                
                tcopy.cdt.insert(candidate);
                tcopy.mark_domain();

                int obtuse_count = tcopy.count_obtuse_triangles();
                
                // If a candidate point eliminates obtuse triangles, store it
                if (obtuse_count < min_obtuse_triangles) {
                    // Update values
                    best_steiner_point = candidate;
                    min_obtuse_triangles = obtuse_count;
                    improved = true;
                }
            }

            if (improved) {
                triangulation.cdt.insert(best_steiner_point);
                triangulation.polygon.push_back(best_steiner_point);
                triangulation.mark_domain();
            }
            // Update the threshold
            triangulation.min_obtuse_triangles = min_obtuse_triangles;
        }
        iterations++;
    }
    
    // Visualize results
    CGAL::draw(triangulation.cdt, triangulation.in_domain);

    // Output final results
    int num_obtuse = triangulation.count_obtuse_triangles();
    int num_steiner = triangulation.cdt.number_of_vertices() - input.num_points;
    std::cout << "Final results: Obtuse triangles: " << num_obtuse
              << ", Steiner points: " << num_steiner << std::endl;

    return triangulation;
}

Triangulation simulated_annealing(const InputJSON& input) {
    // Get parameters from input
    double alpha = input.alpha; // Weight for obtuse triangles in the energy function
    double beta = input.beta;   // Weight for Steiner points in the energy function
    int L = input.L;            // Number of iterations (controls temperature decrement)
    
    // Set initial temperature T <-- 1
    double temperature = 1.0;

    // Initialize triangulation
    Triangulation triangulation = delaunay_const_triangulation(input);

    // Initial energy
    // Energy function combines penalties for obtuse triangles and Steiner points
    int num_obtuse = triangulation.count_obtuse_triangles();
    int num_steiner = triangulation.cdt.number_of_vertices() - input.num_points;
    double energy = alpha * num_obtuse + beta * num_steiner;

    // Random number generator for probabilistic acceptance of "bad" moves
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    // Rate of temperature decrease
    double T_decrement = 1.0 / L; // Decrease temperature by 1/L each iteration

    while (temperature >= 0) {
        // Collecting obtuse triangles
        std::vector<Face_handle> obtuse_faces;
        // Check if the triangle is within the domain and is obtuse
        for (auto face_it = triangulation.cdt.finite_faces_begin(); face_it != triangulation.cdt.finite_faces_end(); ++face_it) {
            // if the triangle is in the domain and is obtuse
            if (triangulation.is_face_in_domain(face_it) && is_obtuse(face_it->vertex(0)->point(), face_it->vertex(1)->point(), face_it->vertex(2)->point())) {
                obtuse_faces.push_back(face_it);
            }
        }

        // If no obtuse triangles remain, stop the optimization
        if (obtuse_faces.empty()) {
            std::cout << "No obtuse triangles left. Optimization complete." << std::endl;
            break;
        }

        // Process each obtuse triangle to attempt improvement
        for (Face_handle fh : obtuse_faces) {
            // Generate candidate Steiner points (5 options)
            // Identify obtuse angle vertex and edges
            Point p0 = fh->vertex(0)->point();
            Point p1 = fh->vertex(1)->point();
            Point p2 = fh->vertex(2)->point();

            // Determine obtuse angle vertex
            Point obtuse_vertex, edge_vertex_1, edge_vertex_2;
            if (angle(p1, p0, p2) == CGAL::OBTUSE) {
                obtuse_vertex = p0;
                edge_vertex_1 = p1;
                edge_vertex_2 = p2;
            } else if (angle(p0, p1, p2) == CGAL::OBTUSE) {
                obtuse_vertex = p1;
                edge_vertex_1 = p0;
                edge_vertex_2 = p2;
            } else {
                obtuse_vertex = p2;
                edge_vertex_1 = p0;
                edge_vertex_2 = p1;
            }

            // Generate candidate points
            std::vector<Point> candidate_points;
            Triangle_2 triangle(p0, p1, p2);
            Point circumcenter = CGAL::circumcenter(triangle); // Circumcenter
            Point centroid = CGAL::centroid(triangle);         // Kentro Varous
            Point projection = project_point_onto_line(obtuse_vertex, edge_vertex_1, edge_vertex_2); // Projection
            Point midpoint = get_midpoint(p0, p1, p2);         // Midpoint

            // Validation checks to stay inbound
            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                && !triangulation.cdt.is_infinite(triangulation.cdt.locate(midpoint))) {
                candidate_points.push_back(midpoint);
            }
            
            if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(circumcenter))) {
                candidate_points.push_back(circumcenter);
            }
            
            if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(centroid))) {
                candidate_points.push_back(centroid);
            }

            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                && !triangulation.cdt.is_infinite(triangulation.cdt.locate(projection))) {
                candidate_points.push_back(projection);
            }

            // Randomly select one candidate point to insert
            std::uniform_int_distribution<size_t> candidate_dist(0, candidate_points.size() - 1);
            Point& candidate = candidate_points[candidate_dist(rng)];

            // Create a copy of the triangulation
            Triangulation tcopy(triangulation);
            tcopy.cdt.insert(candidate);
            tcopy.mark_domain();

            // Calculate new energy
            int new_num_obtuse = tcopy.count_obtuse_triangles();
            int new_num_steiner = tcopy.cdt.number_of_vertices() - input.num_points;
            double new_energy = alpha * new_num_obtuse + beta * new_num_steiner;

            // Calculate the energy difference ΔE
            double delta_energy = new_energy - energy;

            // Metropolis criterion: accept new state probabilistically if ΔE > 0
            if (delta_energy < 0) {
                // Accept
                triangulation.cdt.insert(candidate);
                triangulation.mark_domain();
                triangulation.polygon.push_back(candidate);
                num_obtuse = new_num_obtuse;
                num_steiner = new_num_steiner;
                energy = new_energy;
            } else {
                // Accept with probability e^{-ΔE / T}
                double acceptance_probability = std::exp(-delta_energy / temperature);
                // if acceptance_probability > 1, always accept
                if (uni_dist(rng) <= acceptance_probability) {
                    // Accept worse solution
                    triangulation.cdt.insert(candidate);
                    triangulation.mark_domain();
                    triangulation.polygon.push_back(candidate);
                    num_obtuse = new_num_obtuse;
                    num_steiner = new_num_steiner;
                    energy = new_energy;
                }
            }
        }

        // Decrease temperature
        temperature -= T_decrement;
    }

    // Final output
    std::cout << "Final energy: " << energy << std::endl;
    std::cout << "Obtuse triangles: " << triangulation.count_obtuse_triangles() << ", Steiner points: " << num_steiner << std::endl;

    // Visualize results
    CGAL::draw(triangulation.cdt, triangulation.in_domain);

    return triangulation;
}

Triangulation ant_colony_optimization(const InputJSON& input) {
    // Step 1: Extract parameters from the input JSON
    double alpha = input.alpha;   // Weight for obtuse triangles in energy calculation
    double beta = input.beta;     // Weight for Steiner points in energy calculation
    double xi = input.xi;         // Influence of pheromone in probabilistic selection
    double psi = input.psi;       // Influence of heuristic in probabilistic selection
    double lambda = input.lambda; // Pheromone evaporation rate
    int K = input.kappa;          // Number of ants (colony size)
    int L = input.L;              // Number of optimization cycles

    // Step 2: Define Steiner point types and initialize pheromone values
    enum SteinerPointType {
        VERTEX_PROJECTION = 0, // Projection of the obtuse vertex onto the opposite edge
        CIRCUMCENTER = 1,      // Circumcenter of the obtuse triangle
        MIDPOINT = 2,          // Midpoint of the longest edge of the obtuse triangle
        MEAN_OF_ADJACENT = 3,  // Mean of circumcenters of neighbour obtuse triangles
        NUM_TYPES = 4          // Number of Steiner point types
    };

    // Initialize pheromone values for all Steiner point types (τ_sp > 0 initially)
    std::vector<double> pheromones(NUM_TYPES, 1.0); // τ_sp for each Steiner point type

    // Step 3: Initialize the best triangulation
    Triangulation best_triangulation = delaunay_const_triangulation(input);
    int best_num_obtuse = best_triangulation.count_obtuse_triangles();
    int best_num_steiner = best_triangulation.cdt.number_of_vertices() - input.num_points;
    double best_energy = alpha * best_num_obtuse + beta * best_num_steiner;

    // Step 4: Begin optimization cycles (outer loop: cycles c = 1 to L)
    for (int c = 1; c <= L; ++c) {
        // Each cycle simulates the actions of K ants, who independently improve the triangulation
        std::vector<Triangulation> ant_triangulations(K); // Stores triangulations created by each ant
        std::vector<double> ant_energies(K);              // Stores energy of each ant's triangulation
        std::vector<int> ant_used_options(K, -1);         // Tracks the Steiner point option used by each ant
        
        // Step 5: For each ant k, improve the triangulation by addressing obtuse triangles
        for (int k = 0; k < K; ++k) {
            // **CHANGED HERE**  
            // Instead of a direct assignment (=), do a proper copy so constraints are not lost.
            Triangulation ant_triangulation;
            ant_triangulation.cdt.copy_triangulation(best_triangulation.cdt);
            ant_triangulation.in_domain           = best_triangulation.in_domain;
            ant_triangulation.min_obtuse_triangles = best_triangulation.min_obtuse_triangles;
            ant_triangulation.polygon             = best_triangulation.polygon;
            ant_triangulation.mark_domain();

            // Identify obtuse triangles
            std::vector<Face_handle> obtuse_faces; // List of obtuse triangles
            // Check if the triangle is within the domain and is obtuse
            for (auto face_it = ant_triangulation.cdt.finite_faces_begin(); face_it != ant_triangulation.cdt.finite_faces_end(); ++face_it) {
                if (ant_triangulation.is_face_in_domain(face_it) &&
                    is_obtuse(face_it->vertex(0)->point(), face_it->vertex(1)->point(), face_it->vertex(2)->point())) {
                    obtuse_faces.push_back(face_it);
                }
            }

            // If no obtuse triangles remain, the triangulation is optimal for this ant
            if (obtuse_faces.empty()) {
                // No obtuse triangles, nothing to do
                ant_triangulations[k] = ant_triangulation;
                ant_energies[k] = alpha * 0 + beta * (ant_triangulation.cdt.number_of_vertices() - input.num_points);
                continue;
            }

            // Randomly select one obtuse triangle for improvement
            std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count() + k); // Seed with ant-specific value
            std::uniform_int_distribution<size_t> face_dist(0, obtuse_faces.size() - 1); // Randomly select obtuse triangle
            Face_handle fh = obtuse_faces[face_dist(rng)]; // Selected obtuse triangle

            // For that triangle, compute possible Steiner point options
            Point p0 = fh->vertex(0)->point();
            Point p1 = fh->vertex(1)->point();
            Point p2 = fh->vertex(2)->point();

            // Determine which vertex is the obtuse angle
            Point obtuse_vertex, edge_vertex_1, edge_vertex_2;
            if (CGAL::angle(p1, p0, p2) == CGAL::OBTUSE) {
                obtuse_vertex = p0;
                edge_vertex_1 = p1;
                edge_vertex_2 = p2;
            } else if (CGAL::angle(p0, p1, p2) == CGAL::OBTUSE) {
                obtuse_vertex = p1;
                edge_vertex_1 = p0;
                edge_vertex_2 = p2;
            } else {
                obtuse_vertex = p2;
                edge_vertex_1 = p0;
                edge_vertex_2 = p1;
            }

            // Compute Steiner point candidates and evaluate their heuristic values η.
            // Probabilities P_sp are computed using pheromone values τ and heuristics η.
            // Ant selects a Steiner point probabilistically based on P_sp.
            // Ant inserts the selected Steiner point and evaluates the new triangulation's energy.
            // Ant updates the pheromone values τ based on the energy of the new triangulation.
            std::vector<Point> candidate_points(NUM_TYPES); // Candidate Steiner points for the obtuse triangle

            Triangle_2 triangle(p0, p1, p2); // Triangle formed by the obtuse vertices
            Point circumcenter = CGAL::circumcenter(triangle);
            Point projection = project_point_onto_line(obtuse_vertex, edge_vertex_1, edge_vertex_2);
            Point midpoint = CGAL::midpoint(edge_vertex_1, edge_vertex_2);
            // For MEAN_OF_ADJACENT, we compute the mean of adjacent obtuse triangles' circumcenters
            Point mean_adjacent;
            int adjacent_count = 0;

            // Collect circumcenters of adjacent obtuse triangles
            std::vector<Point> adjacent_circumcenters;
            for (int i = 0; i < 3; ++i) { // Iterate over the three edges of the obtuse triangle
                Face_handle neighbor = fh->neighbor(i); // Get the neighboring triangle
                if (!ant_triangulation.cdt.is_infinite(neighbor) && ant_triangulation.is_face_in_domain(neighbor)) {
                    Point np0 = neighbor->vertex(0)->point();
                    Point np1 = neighbor->vertex(1)->point();
                    Point np2 = neighbor->vertex(2)->point();
                    if (is_obtuse(np0, np1, np2)) { // Check if the neighbor is obtuse
                        Triangle_2 ntriangle(np0, np1, np2);
                        Point ncc = CGAL::circumcenter(ntriangle);
                        adjacent_circumcenters.push_back(ncc);
                        adjacent_count++;
                    }
                }
            }
            // Compute mean of circumcenters if there are at least 2 adjacent obtuse triangles
            if (adjacent_count >= 2) {
                // Compute mean of circumcenters
                K::FT x_sum = 0, y_sum = 0;
                for (const auto& pt : adjacent_circumcenters) {
                    x_sum += pt.x();
                    y_sum += pt.y();
                }
                mean_adjacent = Point(x_sum / adjacent_count, y_sum / adjacent_count); // Mean of circumcenters
            } else {
                // If fewer than 2 adjacent obtuse triangles, we can set mean_adjacent to some default or skip this option
                mean_adjacent = Point(0, 0); // Placeholder
            }

            candidate_points[VERTEX_PROJECTION] = projection;
            candidate_points[CIRCUMCENTER] = circumcenter;
            candidate_points[MIDPOINT] = midpoint;
            candidate_points[MEAN_OF_ADJACENT] = mean_adjacent;

            // Compute radius-to-height ratio ρ
            double circumradius = std::sqrt(CGAL::to_double(CGAL::squared_radius(p0, p1, p2)));

            // Edges
            double lAB = std::sqrt(CGAL::to_double(CGAL::squared_distance(p0, p1)));
            double lBC = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));
            double lCA = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, p0)));

            // Find the longest side length (lmax)
            K::FT lmax = std::max({lAB, lBC, lCA});

            // Compute the area of the triangle
            K::FT area = CGAL::abs(triangle.area());

            // Compute the height relative to the longest side
            K::FT height = (2 * area) / lmax;

            // Avoid division by zero
            if (height == 0) {
                height = 1e-8;
            }

            // Compute ρ
            double rho = CGAL::to_double(circumradius / height); // Radius-to-height ratio

            // Compute η_sp for each option
            std::vector<double> eta(NUM_TYPES, 0.0); // Heuristic values for each Steiner point type

            // VERTEX_PROJECTION
            if (rho > 1.0) {
                eta[VERTEX_PROJECTION] = std::max(0.0, (rho - 1.0) / rho); // Heuristic value for vertex projection
            } else {
                eta[VERTEX_PROJECTION] = 0.0;
            }

            // CIRCUMCENTER
            eta[CIRCUMCENTER] = rho / (2.0 + rho);

            // MIDPOINT
            if (rho < 1.5) {
                eta[MIDPOINT] = std::max(0.0, (3.0 - 2.0 * rho) / 3.0);
            } else {
                eta[MIDPOINT] = 0.0;
            }

            // MEAN_OF_ADJACENT
            if (adjacent_count >= 2) {
                eta[MEAN_OF_ADJACENT] = 1.0;
            } else {
                eta[MEAN_OF_ADJACENT] = 0.0;
            }

            // Now, compute P_sp(k) for each option
            std::vector<double> probabilities(NUM_TYPES, 0.0);
            double sum_probabilities = 0.0;

            for (int sp = 0; sp < NUM_TYPES; ++sp) {
                double tau = pheromones[sp];
                double eta_sp = eta[sp];
                double value = std::pow(tau, xi) * std::pow(eta_sp, psi);
                probabilities[sp] = value;
                sum_probabilities += value;
            }

            // Normalize probabilities
            if (sum_probabilities > 0) {
                for (int sp = 0; sp < NUM_TYPES; ++sp) {
                    probabilities[sp] /= sum_probabilities;
                }
            } else {
                // If sum_probabilities is zero, assign equal probabilities to all options
                for (int sp = 0; sp < NUM_TYPES; ++sp) {
                    probabilities[sp] = 1.0 / NUM_TYPES;
                }
            }

            // Select an option according to probabilities
            std::discrete_distribution<int> option_dist(probabilities.begin(), probabilities.end());
            int selected_option_index = option_dist(rng); // Index of the selected option

            Point selected_point = candidate_points[selected_option_index]; // Selected Steiner point

            // Insert the point into the triangulation
            ant_triangulation.cdt.insert(selected_point);
            ant_triangulation.mark_domain(); // Update domain markers

            // Evaluate the triangulation (compute energy)
            int num_obtuse = ant_triangulation.count_obtuse_triangles();
            int num_steiner = ant_triangulation.cdt.number_of_vertices() - input.num_points;
            double energy = alpha * num_obtuse + beta * num_steiner;

            ant_triangulations[k] = ant_triangulation;
            ant_energies[k] = energy;
            ant_used_options[k] = selected_option_index; // Record the used option
        }

        // Step 6: Update the global best triangulation based on the ants' results.
        for (int k = 0; k < K; ++k) {
            if (ant_energies[k] < best_energy) {
                best_energy = ant_energies[k];
                best_triangulation = ant_triangulations[k];
                best_num_obtuse = best_triangulation.count_obtuse_triangles();
                best_num_steiner = best_triangulation.cdt.number_of_vertices() - input.num_points;
            }
        }

        // Step 7: Update pheromone trails based on the ants' performance.
        // Pheromones are updated based on the number of obtuse triangles reduced by each ant.
        // For each Steiner point option sp (type)
        std::vector<double> delta_pheromones(NUM_TYPES, 0.0);
        for (int k = 0; k < K; ++k) {   // For each ant
            // For the option used by ant k
            int sp = ant_used_options[k];
            if (sp == -1) continue; // No option used
            Triangulation& ant_triangulation = ant_triangulations[k];
            int num_obtuse = ant_triangulation.count_obtuse_triangles();
            int num_steiner = ant_triangulation.cdt.number_of_vertices() - input.num_points;
            if (num_obtuse < best_num_obtuse) {
                // The ant reduced the number of obtuse triangles
                double delta_tau = 1.0 / (1.0 + alpha * num_obtuse + beta * num_steiner); // Δτ_sp, based on energy
                delta_pheromones[sp] += delta_tau; // Accumulate Δτ_sp
            }
        }

        // Update τ_sp
        // Apply pheromone evaporation and add the accumulated Δτ_sp
        for (int sp = 0; sp < NUM_TYPES; ++sp) {
            pheromones[sp] = (1.0 - lambda) * pheromones[sp] + delta_pheromones[sp];
        }
    }

    // Final output
    std::cout << "Final energy: " << best_energy << std::endl;
    std::cout << "Obtuse triangles: " << best_num_obtuse << ", Steiner points: " << best_num_steiner << std::endl;

    // Visualize results
    CGAL::draw(best_triangulation.cdt, best_triangulation.in_domain);
    return best_triangulation;
}
