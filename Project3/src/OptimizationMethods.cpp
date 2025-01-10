#include "../includes/OptimizationMethods.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include <random>
#include <chrono>
#include <cmath>

#include "../includes/draw.h"

Triangulation local_search(const InputJSON& input) {
    std::cout << "Starting local search...\n";
    Triangulation triangulation = delaunay_const_triangulation(input);
    
    bool improved = true;
    int iterations = 0;
    int L = input.L;

    // For detecting deadlock
    int stepsWithoutImprovement = 0;
    
    // For p(n) tracking
    std::vector<double> p_values; 
    p_values.reserve(L + 1);

    int obtuse_prev = triangulation.count_obtuse_triangles();
    int sp_prev     = triangulation.cdt.number_of_vertices();

    // Extra to avoid infinite randomizing
    int loopNoProgressLimit = 2 * input.random_deadlock_threshold; 

    while (improved && (iterations < L)) {
        //Mark facets that are inside the domain bounded by the polygon
        triangulation.mark_domain();

        // Flag to indicate if the triangulation was improved
        bool local_improved = false;

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

            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "midpoint") 
            != input.steiner_methods.end()) {
                // Validation checks to stay inbound
                if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                    && !triangulation.cdt.is_infinite(triangulation.cdt.locate(midpoint))) {
                    std::cout << "Inserting steiner to midpoint" << std::endl;
                    candidate_points.push_back(midpoint);
                }
            }

            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "circumcenter") 
            != input.steiner_methods.end()) {    
                if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(circumcenter))) {
                    std::cout << "Inserting steiner to circumcenter" << std::endl;
                    candidate_points.push_back(circumcenter);
                }
            }
            
            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "centroid") 
            != input.steiner_methods.end()) {    
                if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(centroid))) {
                    std::cout << "Inserting steiner to centroid" << std::endl;
                    candidate_points.push_back(centroid);
                }
            }

            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "projection")
                != input.steiner_methods.end()) {
                std::cout << "projection" << std::endl;
                if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                    && !triangulation.cdt.is_infinite(triangulation.cdt.locate(projection))) {
                    std::cout << "Inserting steiner to projection" << std::endl;
                    candidate_points.push_back(projection);
                }
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
                    local_improved = true;
                }
            }

            if (local_improved) {
                triangulation.cdt.insert(best_steiner_point);
                triangulation.polygon.push_back(best_steiner_point);
                triangulation.mark_domain();
                // Update the threshold
                triangulation.min_obtuse_triangles = min_obtuse_triangles;
                break;
            }
        }

        if (local_improved) {
            stepsWithoutImprovement = 0; // reset
        } else {
            stepsWithoutImprovement += 1;
        }

        improved = local_improved;
        iterations++;

        // Possibly randomize if stuck
        randomize_if_stuck(triangulation, stepsWithoutImprovement, input.random_deadlock_threshold, input.randomization_enabled);
        // If we remain stuck for too long, break (avoid infinite loop)
        if (stepsWithoutImprovement >= loopNoProgressLimit) {
            std::cerr << "[local_search] Breaking out: too many steps without improvement.\n";
            break;
        }

        // measure p(n)
        int obtuse_curr = triangulation.count_obtuse_triangles();
        int sp_curr     = triangulation.cdt.number_of_vertices();
        double pvalue   = compute_pn(obtuse_prev, obtuse_curr, sp_prev, sp_curr);
        p_values.push_back(pvalue);

        obtuse_prev = obtuse_curr;
        sp_prev     = sp_curr;
    }

    // Compute final energy using the same alpha/beta as other methods
    int num_obtuse   = triangulation.count_obtuse_triangles();
    int num_steiner  = triangulation.cdt.number_of_vertices() - input.num_points;
    double final_energy = input.alpha * num_obtuse + input.beta * num_steiner;
    // Visualize results
    // CGAL::draw(triangulation.cdt, triangulation.in_domain);

    // Output final results
    std::cout << "[local_search] Final results:\n"
              << "    Obtuse triangles: " << num_obtuse << "\n"
              << "    Steiner points:   " << num_steiner << "\n"
              << "    Energy:           " << final_energy << std::endl;

    return triangulation;
}

Triangulation simulated_annealing(const InputJSON& input) {
    std::cout << "Starting simulated annealing...\n";
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
    
    // For p(n)
    std::vector<double> p_values;
    p_values.reserve(L + 1);

    // For deadlock detection
    int stepsWithoutImprovement = 0;
    int obtuse_prev = num_obtuse;
    int sp_prev     = triangulation.cdt.number_of_vertices();

    // Extra to avoid infinite randomizing
    int loopNoProgressLimit = 2 * input.random_deadlock_threshold;

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
            std::cout << "[simulated_annealing] No obtuse triangles left. Optimization complete." << std::endl;
            break;
        }

        // Process each obtuse triangle to attempt improvement
        bool improved = false;
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

            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "midpoint") 
            != input.steiner_methods.end()) {
                // Validation checks to stay inbound
                if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                    && !triangulation.cdt.is_infinite(triangulation.cdt.locate(midpoint))) {
                    candidate_points.push_back(midpoint);
                }
            }

            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "circumcenter") 
            != input.steiner_methods.end()) {
                if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(circumcenter))) {
                    candidate_points.push_back(circumcenter);
                }
            }

            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "centroid") 
            != input.steiner_methods.end()) {
                if (!triangulation.cdt.is_infinite(triangulation.cdt.locate(centroid))) {
                    candidate_points.push_back(centroid);
                }
            }

            if (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "projection") 
            != input.steiner_methods.end()) {
                if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)
                    && !triangulation.cdt.is_infinite(triangulation.cdt.locate(projection))) {
                    candidate_points.push_back(projection);
                }
            }
            if (candidate_points.empty()) {
                continue; 
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
                improved = true;
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
                    improved = true;
                }
            }
            // If improved, break so we recalc obtuse set
            if (improved) break;
        }

        if (improved) {
            stepsWithoutImprovement = 0;
        } else {
            stepsWithoutImprovement++;
        }

        // Possibly randomize
        randomize_if_stuck(triangulation, stepsWithoutImprovement, input.random_deadlock_threshold, input.randomization_enabled);

        // If we remain stuck for too long, break
        if (stepsWithoutImprovement >= loopNoProgressLimit) {
            std::cerr << "[simulated_annealing] Breaking out: too many steps w/o improvement.\n";
            break;
        }

        // measure p(n)
        int obtuse_curr = triangulation.count_obtuse_triangles();
        int sp_curr     = triangulation.cdt.number_of_vertices();
        double pvalue   = compute_pn(obtuse_prev, obtuse_curr, sp_prev, sp_curr);
        p_values.push_back(pvalue);

        obtuse_prev = obtuse_curr;
        sp_prev     = sp_curr;

        // Decrease temperature
        temperature -= T_decrement;
        if (temperature < 0) temperature = 0;
    }

    // Compute average p-bar
    double sum_p = 0.0;
    for (double v : p_values) {
        sum_p += v;
    }
    double p_bar = (p_values.empty()) ? 0.0 : sum_p / p_values.size();
    std::cout << "[simulated_annealing] p_bar = " << p_bar << std::endl;
    // Final output
    std::cout << "[simulated_annealing] Final energy: " << energy << std::endl;
    std::cout << "with parameters alpha = " << alpha << " and beta = " << beta << std::endl;
    std::cout << "[simulated_annealing] obtuse: " << triangulation.count_obtuse_triangles()
              << ", steiner: " << (triangulation.cdt.number_of_vertices() - input.num_points) 
              << std::endl;
    // Visualize results
    // CGAL::draw(triangulation.cdt, triangulation.in_domain);

    return triangulation;
}

Triangulation ant_colony_optimization(const InputJSON& input) {
    std::cout << "Starting ant colony optimization...\n";
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

    // For p(n) tracking
    std::vector<double> p_values;
    p_values.reserve(L+1);
    int obtuse_prev = best_num_obtuse;
    int sp_prev     = best_triangulation.cdt.number_of_vertices();

    // track how many cycles we fail to improve
    int cyclesNoImprove = 0;
    int maxNoImprove    = input.random_deadlock_threshold;

    // Step 4: Begin optimization cycles (outer loop: cycles c = 1 to L)
    for (int c = 1; c <= L; ++c) {
        // Each cycle simulates the actions of K ants, who independently improve the triangulation
        std::vector<Triangulation> ant_triangulations(K); // Stores triangulations created by each ant
        std::vector<double> ant_energies(K);              // Stores energy of each ant's triangulation
        std::vector<int> ant_used_options(K, -1);         // Tracks the Steiner point option used by each ant
        // We'll track if any ant improved the best solution
        bool cycleImproved = false;
        // Step 5: For each ant k, improve the triangulation by addressing obtuse triangles
        for (int k = 0; k < K; ++k) {
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

            bool allow_projection   = (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "projection")   != input.steiner_methods.end());
            bool allow_circumcenter = (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "circumcenter") != input.steiner_methods.end());
            bool allow_midpoint     = (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "midpoint")     != input.steiner_methods.end());
            bool allow_meanadj      = (std::find(input.steiner_methods.begin(), input.steiner_methods.end(), "mean_of_adjacent") != input.steiner_methods.end());

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
                // If the user did not allow the corresponding method, zero out its contribution
                if ((sp == VERTEX_PROJECTION && !allow_projection) ||
                    (sp == CIRCUMCENTER       && !allow_circumcenter) ||
                    (sp == MIDPOINT           && !allow_midpoint) ||
                    (sp == MEAN_OF_ADJACENT   && !allow_meanadj))
                {
                    probabilities[sp] = 0.0;
                    continue;
                }
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
                cycleImproved = true;
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

        // measure p(n)
        int obtuse_curr = best_triangulation.count_obtuse_triangles();
        int sp_curr     = best_triangulation.cdt.number_of_vertices();
        double pvalue   = compute_pn(obtuse_prev, obtuse_curr, sp_prev, sp_curr);
        p_values.push_back(pvalue);
        obtuse_prev = obtuse_curr;
        sp_prev     = sp_curr;

        // If no improvement, increment cyclesNoImprove
        if (!cycleImproved) {
            cyclesNoImprove++;
        } else {
            cyclesNoImprove = 0;
        }

        // Possibly do random insertion if stuck:
        // We can treat "cyclesNoImprove" similarly to stepsWithoutImprovement.
        randomize_if_stuck(best_triangulation, cyclesNoImprove, input.random_deadlock_threshold, input.randomization_enabled);

        // If still no improvement for too long, break
        if (cyclesNoImprove >= input.random_deadlock_threshold) {
            std::cerr << "[ant_colony] Breaking out: too many cycles w/o improvement.\n";
            break;
        }
    }

    // Final output
    double sum_p = 0.0;
    for (double v : p_values) sum_p += v;
    double p_bar = (p_values.empty()) ? 0.0 : sum_p / p_values.size();
    std::cout << "[ant_colony] p_bar = " << p_bar << std::endl;
    std::cout << "[ant_colony] final energy: " << best_energy << std::endl;
    std::cout << "[ant_colony] obtuse: " << best_num_obtuse
              << ", steiner: " << best_num_steiner << std::endl;

    // Visualize results
    // CGAL::draw(best_triangulation.cdt, best_triangulation.in_domain);
    return best_triangulation;
}

AutoMethodResult auto_method(const InputJSON& original_input) 
{
    // 1) Detect which category the input belongs to (if you still want to track it)
    std::string category = detect_category(original_input);
    std::cout << "[auto_method] Detected category: " << category << std::endl;

    // 2) We will try ALL methods with multiple param combinations, 
    //    ignoring the category-based selection. 
    //    The category can, if you wish, prune some parameter sets or methods, 
    //    but the assignment says we must check all. So let's do it.

    Triangulation best_triang;
    InputJSON     best_input;  // Will hold the method/params for the best triang
    double        best_energy  = std::numeric_limits<double>::infinity();
    bool          found_best   = false;

    // Helper lambda to evaluate energy of a triangulation
    auto compute_energy = [&](const Triangulation& t, const InputJSON& in) {
        int obtuse  = t.count_obtuse_triangles();
        int steiner = t.cdt.number_of_vertices() - in.num_points;
        return (in.alpha * obtuse) + (in.beta * steiner);
    };

    // Make local copies of input to override parameters on each test
    // so we don't mutate the original.
    auto run_local_search = [&](int I) {
        InputJSON tmp = original_input;
        tmp.method = "local";
        // We store the local search iteration limit in tmp.L or in parameters map
        tmp.parameters["I"] = boost::json::value(I);
        tmp.L = I;  // If your local_search code reads from L
        tmp.parameters["L"] = I;
        Triangulation triang = local_search(tmp);

        double energy = compute_energy(triang, tmp);
        if (energy < best_energy) {
            best_energy    = energy;
            best_triang    = triang;
            best_input   = tmp;
            found_best     = true;
        }
    };

    auto run_sim_anneal = [&](double a, double b, int Lval) {
        InputJSON tmp = original_input;
        tmp.method = "sa";
        tmp.alpha  = a;
        tmp.parameters["alpha"] = a;
        tmp.beta   = b;
        tmp.parameters["beta"] = b;
        tmp.L      = Lval;
        tmp.parameters["L"] = Lval;
        Triangulation triang = simulated_annealing(tmp);

        double energy = compute_energy(triang, tmp);
        if (energy < best_energy) {
            best_energy = energy;
            best_triang = triang;
            best_input   = tmp;
            found_best  = true;
        }
    };

    auto run_ant_colony = [&](double a, double b, double x, double p, double lam, int kappa_val, int Lval) {
        InputJSON tmp = original_input;
        tmp.method  = "ant";
        tmp.alpha   = a;
        tmp.parameters["alpha"] = a;
        tmp.beta    = b;
        tmp.parameters["beta"] = b;
        tmp.xi      = x;
        tmp.parameters["xi"] = x;
        tmp.psi     = p;
        tmp.parameters["psi"] = p;
        tmp.lambda  = lam;
        tmp.parameters["lambda"] = lam;
        tmp.kappa   = kappa_val;
        tmp.parameters["kappa"] = kappa_val;
        tmp.L       = Lval;
        tmp.parameters["L"] = Lval;

        Triangulation triang = ant_colony_optimization(tmp);
        double energy = compute_energy(triang, tmp);
        if (energy < best_energy) {
            best_energy = energy;
            best_triang = triang;
            best_input   = tmp;
            found_best  = true;
        }
    };

    // ------------------------------------------------------------------------
    // 3) Enumerate all parameter sets for each method
    // ------------------------------------------------------------------------

    // =========== Local Search =============
    {
        std::vector<int> local_I_values = {50, 100, 200};
        for (int I : local_I_values) {
            run_local_search(I);
        }
    }

    // =========== Simulated Annealing =============
    {
        std::vector<std::pair<double,double>> alphaBetaSA = {
            {1.0, 1.0},
            {1.0, 2.0},
            {2.0, 2.0}
        };
        std::vector<int> Lvals = {1000, 5000};

        for (auto &ab : alphaBetaSA) {
            for (auto &Lval : Lvals) {
                run_sim_anneal(ab.first, ab.second, Lval);
            }
        }
    }

    // =========== Ant Colony =============
    {
        std::vector<std::pair<double,double>> alphaBetaAnt = {
            {1.0, 1.0},
            {1.0, 2.0}
        };
        std::vector<std::pair<double,double>> xiPsi = {
            {1.0, 1.0},
            {2.0, 3.0}
        };
        std::vector<double> lambdas = {0.3, 0.5};
        std::vector<int> kappaVals = {5, 10};
        std::vector<int> antLvals  = {30, 50};

        for (auto &ab : alphaBetaAnt) {
            for (auto &xp : xiPsi) {
                for (auto &lam : lambdas) {
                    for (auto &kv : kappaVals) {
                        for (auto &Lv : antLvals) {
                            run_ant_colony(ab.first, ab.second, xp.first, xp.second, lam, kv, Lv);
                        }
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 4) Possibly do a final randomization if too many obtuse remain
    // ------------------------------------------------------------------------
    if (found_best) 
    {
        int final_obtuse = best_triang.count_obtuse_triangles();
        int final_steiner = best_triang.cdt.number_of_vertices() - original_input.num_points;
        double final_energy = original_input.alpha * final_obtuse + original_input.beta * final_steiner;

        // Extract the best method and Steiner options
        std::string best_method = best_input.method; // e.g., "local", "sa", or "ant"
        std::vector<std::string> steiner_options = best_input.steiner_methods; // List of Steiner options used

        std::cout << "[auto_method] Final Optimization State:\n"
                << "    Method:           " << best_method << "\n"
                << "    Steiner Options:  ";
        for (const auto& option : steiner_options) {
            std::cout << option << " ";
        }
        std::cout << "\n"
                << "    Obtuse triangles: " << final_obtuse << "\n"
                << "    Steiner points:   " << final_steiner << "\n"
                << "    Energy:           " << final_energy << "\n";
        // Example heuristic
        if (final_obtuse > 5) {
            std::cout << "[auto_method] Too many obtuse triangles in the final triangulation. Randomizing...\n";
            best_triang.randomizationUsed = true;
            std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
            
            // Rough bounding box or you can compute from region_boundary
            double xmin = 0.0, xmax = 1000.0, ymin = 0.0, ymax = 1000.0;
            // or call compute_bounding_box(original_input, xmin, xmax, ymin, ymax);

            std::uniform_real_distribution<double> distx(xmin, xmax);
            std::uniform_real_distribution<double> disty(ymin, ymax);

            for (int i = 0; i < 5; i++) {
                Point rp(distx(rng), disty(rng));
                best_triang.cdt.insert(rp);
                best_triang.mark_domain();
            }
        }
        best_triang.randomizationUsed = true;
        best_triang.mark_domain();
        AutoMethodResult result;
        result.triang     = best_triang;
        result.best_input = best_input;  // Holds the correct method/params
        return result;
    }
    else {
        // If somehow we found no triang (very unlikely), fallback
        std::cerr << "[auto_method] Warning: no method produced a valid triang. Returning Delaunay.\n";
                Triangulation fallback = delaunay_const_triangulation(original_input);
        AutoMethodResult result;
        result.triang     = fallback;
        result.best_input = original_input;  // Or set method="none" if you prefer
        return result;
    }
}
