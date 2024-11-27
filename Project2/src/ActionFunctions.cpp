#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include <queue>
#include <cmath>

// Function to check if a triangle is obtuse
bool is_obtuse(const Point& A, const Point& B, const Point& C) {

    if (CGAL::angle(A, B, C) == CGAL::OBTUSE || CGAL::angle(B, C, A) == CGAL::OBTUSE || CGAL::angle(C, A, B) == CGAL::OBTUSE) {
        return true;
    }
    return false;
}

// Function to try edge flips
void edge_flip(Triangulation& triangulation) {

    bool flipped = true;
    while (flipped) {
        flipped = false;
        // Iterate through all triangles
        for (CDT::Finite_faces_iterator fit = triangulation.cdt.finite_faces_begin(); fit != triangulation.cdt.finite_faces_end(); ++fit) {
            Face_handle fh = fit;

            // Check if the triangle is obtuse
            if (is_obtuse(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point())) {
                // For each neighbor
                for (int i = 0; i < 3; ++i) {
                    if (triangulation.cdt.is_flipable(fh, i)) {
                        // Count obtuse triangles before the flip
                        int before_flip = triangulation.count_obtuse_triangles();

                        // Flip the edge
                        triangulation.cdt.flip(fh, i);

                        // Count obtuse triangles after the flip
                        int after_flip = triangulation.count_obtuse_triangles();

                        // Flip back if the number of obtuse triangles has not decreased
                        if (after_flip >= before_flip) {
                            triangulation.cdt.flip(fh, i);
                        } else {
                            flipped = true;
                            break;
                        }
                    }
                }
            }
        }
    }
}

// Function to project point A onto the line defined by B and C
Point project_point_onto_line(const Point& A, const Point& B, const Point& C) {

    // Vector from B to C
    K::Vector_2 BC = C - B;
    // Vector from B to A
    K::Vector_2 BA = A - B;

    // Scalar projection of BA onto BC
    K::FT t = (BA * BC) / (BC * BC);

    // Clamp t to the interval [0, 1] to ensure the projection falls on the line segment BC
    t = std::max(K::FT(0), std::min(K::FT(1), t));

    // Compute the projection point using B + t * (C - B)
    return B + t * BC;
}

// Function to check if an edge is inbound
bool is_edge_valid(const Point& v1, const Point& v2, const Polygon_2& polygon) {
    
    // Check if both vertices are inside the polygon or on the boundary
    CGAL::Bounded_side v1_side = polygon.bounded_side(v1);
    CGAL::Bounded_side v2_side = polygon.bounded_side(v2);

    if (v1_side == CGAL::ON_UNBOUNDED_SIDE || v2_side == CGAL::ON_UNBOUNDED_SIDE) {
        // If either vertex is outside the polygon, the edge is invalid
        return false;
    }

    // Edge is valid if vertices are valid
    return true;
}

// Function to find the midpoint of the largest edge formed by the points given
Point get_midpoint(const Point& A, const Point& B, const Point& C) {
    // Calculate squared distances of the edges
    K::FT lAB = CGAL::squared_distance(A, B);
    K::FT lBC = CGAL::squared_distance(B, C);
    K::FT lCA = CGAL::squared_distance(C, A);

    // Calculate max and min edges
    K::FT lmax = std::max({lAB, lBC, lCA});
    K::FT lmin = std::min({lAB, lBC, lCA});

    // Calculate the midpoint of the largest edge
    if (lmax == lAB) {
        return CGAL::midpoint(A, B);
    } else if (lmax == lBC) {
        return CGAL::midpoint(B, C);
    } else {
        return CGAL::midpoint(C, A);
    }
}

// Function to add Steiner points
bool add_optimal_steiner(Triangulation& triangulation) {
    
    //Mark facets that are inside the domain bounded by the polygon
    triangulation.mark_domain();

    // Flag to indicate if the triangulation was improved
    bool improved = false;

    Point best_steiner_point;
    int min_obtuse_triangles = triangulation.min_obtuse_triangles;

    // Store triangluation states
    std::priority_queue<TriangulationState> pq;    
    
    // Data structures to store points
    std::vector<Point> initial_steiner_points;
    initial_steiner_points.assign(triangulation.polygon.begin(), triangulation.polygon.end());
    std::vector<Point> current_steiner_points;

    // Store the initial obtuse triangles
    std::vector<Face_handle> obtuse_faces;
    for (auto face_it = triangulation.cdt.finite_faces_begin(); face_it != triangulation.cdt.finite_faces_end(); ++face_it) {
        if (is_obtuse(face_it->vertex(0)->point(), face_it->vertex(1)->point(), face_it->vertex(2)->point())) {
            obtuse_faces.push_back(face_it);
        }
    }

    // Set the initial triangulation as a base line
    pq.push({initial_steiner_points, min_obtuse_triangles});

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

        std::set<Edge> shared_edges;
        // Ensure given triangle does not have constrained edges
        bool all_edges_unconstrained = true;
        for (int i = 0; i < 3; i++) {
            if(triangulation.cdt.is_constrained(Edge(fh, i))) {
                all_edges_unconstrained = false;
                break;
            }
        }

        // Get neighbors that do not have constrained edges
        std::set<Face_handle> valid_faces;
        if (all_edges_unconstrained) {
            // Store valid triangles
            valid_faces.insert(fh);

            // Iterate through each edge of the triangle to examine neighbors
            for (int i = 0; i < 3; i++) {
                Face_handle neighbor_fh = fh->neighbor(i);

                // Check if the neighbor is finite and obtuse
                bool neighbor_edges_unconstrained = true;
                if (!triangulation.cdt.is_infinite(neighbor_fh) 
                    && is_obtuse(neighbor_fh->vertex(0)->point(), neighbor_fh->vertex(1)->point(), neighbor_fh->vertex(2)->point())) {
                    for (int j = 0; j < 3; j++) {
                        Edge edge(neighbor_fh, j);
                        if (triangulation.cdt.is_constrained(edge)) {
                            neighbor_edges_unconstrained = false;
                            break;
                        }
                    }
                }
                if (neighbor_edges_unconstrained) {
                    valid_faces.insert(neighbor_fh);

                    // Store shared edges
                    for (int j = 0; j < 3; j++) {
                        Vertex_handle v1 = fh->vertex(j);
                        Vertex_handle v2 = fh->vertex((j + 1) % 3);

                        // Check if the neighbor has the same edge
                        for (int k = 0; k < 3; k++) {
                            Vertex_handle nv1 = neighbor_fh->vertex(k);
                            Vertex_handle nv2 = neighbor_fh->vertex((k + 1) % 3);

                            if ((v1 == nv1 && v2 == nv2) || (v1 == nv2 && v2 == nv1)) {
                                shared_edges.insert(Edge(fh, j));
                                break;
                            }
                        }
                    }
                }
            }

            // Store the formed polygon
            Polygon_2 polygon_vertices;
            for (const auto&face : valid_faces) {
                for (int i = 0; i < 3; i++) {
                    Point vertex = face->vertex(i)->point();

                    // Ignore duplicates
                    if (std::find(polygon_vertices.begin(), polygon_vertices.end(), vertex) == polygon_vertices.end()) {
                        polygon_vertices.push_back(vertex);
                    }
                }
            }
            if (polygon_vertices.is_convex()) {
                std::vector<Constraint_id> added_constraints;
                // Constrain the polygon
                for (size_t i = 0; i < polygon_vertices.size(); ++i) {
                    Point v1 = polygon_vertices[i];
                    Point v2 = polygon_vertices[(i+1) % polygon_vertices.size()];

                    // Skip if the edge is in shared_edges
                    bool is_shared = false;
                    for (const auto &edge : shared_edges) {
                        Vertex_handle ev1 = edge.first->vertex(edge.second);
                        Vertex_handle ev2 = edge.first->vertex((edge.second + 1) % 3);

                        if ((ev1->point() == v1 && ev2->point() == v2) || (ev1->point() == v2 && ev2->point() == v1)) {
                            is_shared = true;
                            break;
                        }
                    }

                    if (!is_shared) {
                        Constraint_id cid = triangulation.cdt.insert_constraint(v1, v2); // Add boundary constraint
                        added_constraints.push_back(cid);
                    }
                }
                
                std::vector<Point> shared_edge_points;
                
                for (const auto &edge : shared_edges) {
                    Vertex_handle v1 = edge.first->vertex(edge.second);
                    Vertex_handle v2 = edge.first->vertex((edge.second + 1) % 3);

                    shared_edge_points.push_back(v1->point());
                    shared_edge_points.push_back(v2->point());

                    triangulation.remove_no_flip(v1);
                    triangulation.remove_no_flip(v2);
                }

                polygon_centroid = CGAL::centroid(polygon_vertices.begin(), polygon_vertices.end());
                triangulation.cdt.insert(polygon_centroid);

                // Re-insert the stored points
                for (const auto &point : shared_edge_points) {
                    triangulation.cdt.insert(point);
                }
                
                if (!added_constraints.empty()) {
                    for (size_t i = 0; i < added_constraints.size(); ++i) {
                        triangulation.cdt.remove_constraint(added_constraints[i]);
                    }
                }
                std::cout<<"Used"<<std::endl;
            }
        }

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
            
            tcopy.mark_domain();
            tcopy.cdt.insert(candidate);
            tcopy.mark_domain();

            int obtuse_count = tcopy.count_obtuse_triangles();
            
            // If a candidate point eliminates obtuse triangles, store it
            if (obtuse_count < min_obtuse_triangles) {
                // Update values
                best_steiner_point = candidate;
                min_obtuse_triangles = obtuse_count;
                current_steiner_points.push_back(best_steiner_point);
                pq.push({current_steiner_points, obtuse_count});
                
                improved = true;
            }
        }
    }

    // Retrieve the state with the least obtuse triangles
    TriangulationState best_state = pq.top();

    // If the Steiner point addition improved the triangulation
    if (improved) {
        
        // Create another test copy
        Triangulation tcopy2(triangulation);
        tcopy2.mark_domain();
        // Add the Steiner points
        for (auto i = 0; i < best_state.steiner_points.size(); ++i) {
            tcopy2.cdt.insert(best_state.steiner_points[i]);
            tcopy2.polygon.push_back(best_state.steiner_points[i]);
        }

        // Check if the addition of points affects the number of obtuse triangles afterwards
        if (tcopy2.count_obtuse_triangles() < triangulation.count_obtuse_triangles()) {
            for (auto i = 0; i < best_state.steiner_points.size(); ++i) {
                triangulation.cdt.insert(best_state.steiner_points[i]);
                triangulation.polygon.push_back(best_state.steiner_points[i]);
                triangulation.mark_domain();
            }
            
            // Update the threshold
            triangulation.min_obtuse_triangles = best_state.obtuse_triangle_count;
        }

    }

    return improved;
}

// Function to eliminate obtuse triangles
void eliminate_obtuse_triangles(Triangulation& triangulation, int L) {
    bool improved = true;
    int iterations = 0;

    // Loop terminates when adding steiner points doesn't improve the triangulations
    while (improved && iterations < L) {
        improved = add_optimal_steiner(triangulation);
        iterations++;
    }

    // Result print
    std::cout << "Total obtuse triangles: " << triangulation.count_obtuse_triangles() << std::endl;
}