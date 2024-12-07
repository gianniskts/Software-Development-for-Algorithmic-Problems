#include "../includes/OptimizationMethods.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include <random>
#include <chrono>
#include <cmath>

Triangulation local_search(const InputJSON& input) {
    return delaunay_const_triangulation(input);
}

Triangulation simulated_annealing(const InputJSON& input) {
    // Get parameters from input
    double alpha = input.alpha;
    double beta = input.beta;
    int L = input.L;
    // Set initial temperature T <-- 1
    double temperature = 1.0;

    // Initialize triangulation
    Triangulation triangulation = delaunay_const_triangulation(input);

    // Calculate initial energy
    int num_obtuse = triangulation.count_obtuse_triangles();
    int num_steiner = triangulation.cdt.number_of_vertices() - input.num_points;
    double energy = alpha * num_obtuse + beta * num_steiner;

    // Set up random number generator
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

    double T_decrement = 1.0 / L; // Decrease temperature by 1/L each iteration

    while (temperature > 0) {
        // Collect obtuse triangles
        std::vector<Face_handle> obtuse_faces;
        // Iterate through all triangles
        for (auto face_it = triangulation.cdt.finite_faces_begin(); face_it != triangulation.cdt.finite_faces_end(); ++face_it) {
            // if the triangle is in the domain and is obtuse
            if (triangulation.is_face_in_domain(face_it) && is_obtuse(face_it->vertex(0)->point(), face_it->vertex(1)->point(), face_it->vertex(2)->point())) {
                obtuse_faces.push_back(face_it);
            }
        }

        // If no obtuse triangles, break the loop
        if (obtuse_faces.empty()) {
            std::cout << "No obtuse triangles left. Optimization complete." << std::endl;
            break;
        }

        // For each obtuse triangle
        for (Face_handle fh : obtuse_faces) {
            // Generate candidate Steiner points (5 options)
            Point p0 = fh->vertex(0)->point();
            Point p1 = fh->vertex(1)->point();
            Point p2 = fh->vertex(2)->point();

            // Determine which vertex is the obtuse angle
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
            Point circumcenter = CGAL::circumcenter(triangle);
            Point centroid = CGAL::centroid(triangle);
            Point projection = project_point_onto_line(obtuse_vertex, edge_vertex_1, edge_vertex_2);
            Point midpoint = get_midpoint(p0, p1, p2);

            candidate_points.push_back(circumcenter);
            candidate_points.push_back(centroid);
            candidate_points.push_back(projection);
            candidate_points.push_back(midpoint);

            // Add a random point within the triangle
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            double r1 = dist(rng);
            double r2 = dist(rng);
            if (r1 + r2 > 1) {
                r1 = 1 - r1;
                r2 = 1 - r2;
            }
            Point random_point = p0 + (p1 - p0) * r1 + (p2 - p0) * r2;
            candidate_points.push_back(random_point);

            // Randomly select one candidate
            std::uniform_int_distribution<size_t> candidate_dist(0, candidate_points.size() - 1);
            const Point& candidate = candidate_points[candidate_dist(rng)];

            // Create a copy of the triangulation
            Triangulation tcopy(triangulation);
            tcopy.cdt.insert(candidate);
            tcopy.mark_domain();

            // Calculate new energy
            int new_num_obtuse = tcopy.count_obtuse_triangles();
            int new_num_steiner = tcopy.cdt.number_of_vertices() - input.num_points;
            double new_energy = alpha * new_num_obtuse + beta * new_num_steiner;

            // Calculate ΔE
            double delta_energy = new_energy - energy;

            // Decide whether to accept
            if (delta_energy < 0) {
                // Accept
                triangulation = tcopy;
                num_obtuse = new_num_obtuse;
                num_steiner = new_num_steiner;
                energy = new_energy;
            } else {
                // Accept with probability e^{-ΔE / T}
                double acceptance_probability = std::exp(-delta_energy / temperature);
                // if acceptance_probability > 1, always accept
                if (uni_dist(rng) < acceptance_probability) {
                    // Accept worse solution
                    triangulation = tcopy;
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
    std::cout << "Obtuse triangles: " << num_obtuse << ", Steiner points: " << num_steiner << std::endl;

    return triangulation;
}

Triangulation ant_colony_optimization(const InputJSON& input) {
    // Get parameters from input
    double alpha = input.alpha;   // α
    double beta = input.beta;     // β
    double xi = input.xi;         // χ
    double psi = input.psi;       // ψ
    double lambda = input.lambda; // λ
    int K = input.kappa;          // Number of ants
    int L = input.L;              // Number of cycles

    // Initialize pheromone values τ0 > 0 for all Steiner point options (types)
    enum SteinerPointType {
        VERTEX_PROJECTION = 0,
        CIRCUMCENTER = 1,
        MIDPOINT = 2,
        MEAN_OF_ADJACENT = 3,
        NUM_TYPES = 4
    };

    std::vector<double> pheromones(NUM_TYPES, 1.0); // τ_sp for each Steiner point type

    // Initialize the best triangulation
    Triangulation best_triangulation = delaunay_const_triangulation(input);
    int best_num_obtuse = best_triangulation.count_obtuse_triangles();
    int best_num_steiner = best_triangulation.cdt.number_of_vertices() - input.num_points;
    double best_energy = alpha * best_num_obtuse + beta * best_num_steiner;

    // For cycle c = 1 to L
    for (int c = 1; c <= L; ++c) {
        // We'll store the results for each ant
        // For ant k = 1 to K
        std::vector<Triangulation> ant_triangulations(K);
        std::vector<double> ant_energies(K, std::numeric_limits<double>::infinity());
        std::vector<int> ant_used_options(K, -1); // Keep track of options used by each ant
        std::vector<Point> ant_inserted_points(K, Point(0,0));
        std::vector<int> ant_obtuse_counts(K, 0);
        std::vector<int> ant_steiner_counts(K, 0);

        // For each ant, test a solution on a copy
        for (int k = 0; k < K; ++k) {
            Triangulation test_triangulation = best_triangulation;

            // Identify obtuse triangles
            std::vector<Face_handle> obtuse_faces;
            for (auto face_it = test_triangulation.cdt.finite_faces_begin(); 
                 face_it != test_triangulation.cdt.finite_faces_end(); ++face_it) {
                if (test_triangulation.is_face_in_domain(face_it) &&
                    is_obtuse(face_it->vertex(0)->point(), face_it->vertex(1)->point(), face_it->vertex(2)->point())) {
                    obtuse_faces.push_back(face_it);
                }
            }

            if (obtuse_faces.empty()) {
                // No obtuse triangles, nothing to do
                int num_obtuse = test_triangulation.count_obtuse_triangles();
                int num_steiner = test_triangulation.cdt.number_of_vertices() - input.num_points;
                double energy = alpha * num_obtuse + beta * num_steiner;
                ant_triangulations[k] = test_triangulation;
                ant_energies[k] = energy;
                ant_obtuse_counts[k] = num_obtuse;
                ant_steiner_counts[k] = num_steiner;
                continue;
            }

            // Randomly select one obtuse triangle
            std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count() + k);
            std::uniform_int_distribution<size_t> face_dist(0, obtuse_faces.size() - 1);
            Face_handle fh = obtuse_faces[face_dist(rng)];

            // Extract points of the chosen obtuse triangle
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

            // Compute candidate Steiner points
            Triangle_2 triangle(p0, p1, p2);
            Point circumcenter = CGAL::circumcenter(triangle);
            Point projection = project_point_onto_line(obtuse_vertex, edge_vertex_1, edge_vertex_2);
            Point midpoint = CGAL::midpoint(edge_vertex_1, edge_vertex_2);
            // For MEAN_OF_ADJACENT, we need to compute the mean of adjacent obtuse triangles' circumcenters
            Point mean_adjacent;
            int adjacent_count = 0;
            {
                // Collect circumcenters of adjacent obtuse triangles
                std::vector<Point> adjacent_circumcenters;
                for (int i = 0; i < 3; ++i) {
                    Face_handle neighbor = fh->neighbor(i);
                    if (!test_triangulation.cdt.is_infinite(neighbor) && test_triangulation.is_face_in_domain(neighbor)) {
                        Point np0 = neighbor->vertex(0)->point();
                        Point np1 = neighbor->vertex(1)->point();
                        Point np2 = neighbor->vertex(2)->point();
                        if (is_obtuse(np0, np1, np2)) {
                            Triangle_2 ntriangle(np0, np1, np2);
                            Point ncc = CGAL::circumcenter(ntriangle);
                            adjacent_circumcenters.push_back(ncc);
                            adjacent_count++;
                        }
                    }
                }

                if (adjacent_count >= 2) {
                    // Compute mean of circumcenters
                    K::FT x_sum = 0, y_sum = 0;
                    for (const auto& pt : adjacent_circumcenters) {
                        x_sum += pt.x();
                        y_sum += pt.y();
                    }
                    mean_adjacent = Point(x_sum / adjacent_count, y_sum / adjacent_count);
                } else {
                    // If fewer than 2 adjacent obtuse triangles, we can set mean_adjacent to some default or skip this option
                    mean_adjacent = Point(0, 0); // Default if not enough adjacent
                }
            }

            std::vector<Point> candidate_points(NUM_TYPES);
            candidate_points[VERTEX_PROJECTION] = projection;
            candidate_points[CIRCUMCENTER] = circumcenter;
            candidate_points[MIDPOINT] = midpoint;
            candidate_points[MEAN_OF_ADJACENT] = mean_adjacent;

            // Compute radius-to-height ratio ρ
            double circumradius = std::sqrt(CGAL::to_double(CGAL::squared_radius(p0, p1, p2)));

            // Compute the length of each edge
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
            double rho = CGAL::to_double(circumradius / height);

            // Compute η_sp for each option
            std::vector<double> eta(NUM_TYPES, 0.0);

            // VERTEX_PROJECTION
            if (rho > 1.0) {
                eta[VERTEX_PROJECTION] = std::max(0.0, (rho - 1.0) / rho);
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
            int selected_option_index = option_dist(rng);

            Point selected_point = candidate_points[selected_option_index];

            // Test insertion on the copy (test_triangulation)
            test_triangulation.cdt.insert(selected_point);
            test_triangulation.mark_domain();

            // Evaluate the triangulation (compute energy)
            int num_obtuse = test_triangulation.count_obtuse_triangles();
            int num_steiner = test_triangulation.cdt.number_of_vertices() - input.num_points;
            double energy = alpha * num_obtuse + beta * num_steiner;

            // Store ant results
            ant_triangulations[k] = test_triangulation;
            ant_energies[k] = energy;
            ant_used_options[k] = selected_option_index;
            ant_inserted_points[k] = selected_point;
            ant_obtuse_counts[k] = num_obtuse;
            ant_steiner_counts[k] = num_steiner;
        }

        // SaveBestTriangulation[c]: select the best among ants
        // After all ants: choose the best among them
        double cycle_best_energy = best_energy;
        int cycle_best_ant = -1;
        for (int k = 0; k < K; ++k) {
            if (ant_energies[k] < cycle_best_energy) {
                cycle_best_energy = ant_energies[k];
                cycle_best_ant = k;
            }
        }

        // If there's an improvement, apply it to the best_triangulation
        if (cycle_best_ant != -1 && ant_energies[cycle_best_ant] < best_energy) {
            // Insert the chosen point of that ant into best_triangulation
            Point chosen_point = ant_inserted_points[cycle_best_ant];
            best_triangulation.cdt.insert(chosen_point);
            best_triangulation.mark_domain();

            // Update best_energy and counts
            best_energy = ant_energies[cycle_best_ant];
            best_num_obtuse = best_triangulation.count_obtuse_triangles();
            best_num_steiner = best_triangulation.cdt.number_of_vertices() - input.num_points;
        }

        // Update pheromones
        // We'll reward ants that produced solutions with fewer obtuse triangles
        // than the current best_num_obtuse, as in the original logic.
        std::vector<double> delta_pheromones(NUM_TYPES, 0.0);
        for (int k = 0; k < K; ++k) {
            // For the option used by ant k
            int sp = ant_used_options[k];
            if (sp == -1) continue;

            int num_obtuse = ant_obtuse_counts[k];
            int num_steiner = ant_steiner_counts[k];

            // If the ant reduced the number of obtuse triangles relative to best_num_obtuse
            // at the time of evaluation, or contributed to a better energy, give pheromone:
            if (num_obtuse < best_num_obtuse) {
                // The ant reduced the number of obtuse triangles
                double delta_tau = 1.0 / (1.0 + alpha * num_obtuse + beta * num_steiner);
                delta_pheromones[sp] += delta_tau;
            }
        }

        // Update τ_sp
        for (int sp = 0; sp < NUM_TYPES; ++sp) {
            pheromones[sp] = (1.0 - lambda) * pheromones[sp] + delta_pheromones[sp];
        }
    }

    // Final output
    std::cout << "Final energy: " << best_energy << std::endl;
    std::cout << "Obtuse triangles: " << best_num_obtuse << ", Steiner points: " << best_num_steiner << std::endl;

    return best_triangulation;
}
