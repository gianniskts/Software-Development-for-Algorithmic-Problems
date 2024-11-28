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
    // For now, returning a basic triangulation
    return delaunay_const_triangulation(input);
}
