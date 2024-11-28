#include "../includes/OptimizationMethods.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include <random>
#include <chrono>
#include <cmath>

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
            if (is_obtuse(face_it->vertex(0)->point(), face_it->vertex(1)->point(), face_it->vertex(2)->point())) {
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

            // Validation check to stay inbound
            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)) {
                candidate_points.push_back(midpoint);
            }
            
            // If circumcenter is infinite, then get the centroid
            if (!(triangulation.cdt.is_infinite(triangulation.cdt.locate(circumcenter)))) {
                candidate_points.push_back(circumcenter);
            } else {
                candidate_points.push_back(centroid);
            }

            // If inbound, get the projection
            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)) {
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
    
    return triangulation;
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

            // Validation check to stay inbound
            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)) {
                candidate_points.push_back(midpoint);
            }
            
            // If circumcenter is infinite, then get the centroid
            if (!(triangulation.cdt.is_infinite(triangulation.cdt.locate(circumcenter)))) {
                candidate_points.push_back(circumcenter);
            } else {
                candidate_points.push_back(centroid);
            }

            // If inbound, get the projection
            if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)) {
                candidate_points.push_back(projection);
            }

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
    std::cout << "Obtuse triangles: " << triangulation.count_obtuse_triangles() << ", Steiner points: " << num_steiner << std::endl;

    return triangulation;
}

Triangulation ant_colony_optimization(const InputJSON& input) {
    // For now, returning a basic triangulation
    return delaunay_const_triangulation(input);
}
