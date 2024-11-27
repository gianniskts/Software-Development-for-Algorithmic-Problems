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
    bool polygon_method_used = false;

    Point best_steiner_point;
    int min_obtuse_triangles = triangulation.min_obtuse_triangles;   
    
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
                }
            }

            // Store the formed polygon
            Polygon_2 polygon_points;
            std::vector<Vertex_handle> polygon_vertices;
            for (const auto&face : valid_faces) {
                for (int i = 0; i < 3; i++) {
                    Point vertex = face->vertex(i)->point();

                    // Skip duplicates
                    if (std::find(polygon_points.begin(), polygon_points.end(), vertex) == polygon_points.end()) {
                        polygon_points.push_back(vertex);
                        polygon_vertices.push_back(face->vertex(i));
                    }
                }
            }
            // Check if the polygon formed is convex
            if (polygon_points.is_convex()) {
                std::set<Constraint_id> added_constraints;
                Triangulation poly_copy(triangulation);
                // Constrain the polygon
                for (auto edge = polygon_points.edges_begin(); edge != polygon_points.edges_end(); ++edge) {
                    Segment_2 current_edge(*edge);
                    if (polygon_points.bounded_side(current_edge.source()) == CGAL::ON_BOUNDARY &&
                        polygon_points.bounded_side(current_edge.target()) == CGAL::ON_BOUNDARY) {
                        
                        // Store cid to delete it later
                        Constraint_id cid = triangulation.cdt.insert_constraint(current_edge.source(), current_edge.target());
                        added_constraints.insert(cid);
                    }
                }
                
                // Remove the vertices of the polygon from the triangulation
                for (size_t i = 0 ; i < polygon_vertices.size(); ++i) {
                    Vertex_handle vh = polygon_vertices[i];
                    // Make sure a boundary point is not removed
                    if (!triangulation.are_there_incident_constraints(vh)) {
                        triangulation.remove_no_flip(vh);
                    }
                }

                // Add the steiner point at the center of the polygon
                polygon_centroid = CGAL::centroid(polygon_points.begin(), polygon_points.end());
                triangulation.cdt.insert(polygon_centroid);
                triangulation.mark_domain();
                triangulation.polygon.push_back(polygon_centroid);
                
                
                // Re-insert the stored points
                for (const auto &point : polygon_points) {
                    triangulation.cdt.insert(point);
                    triangulation.mark_domain();
                }
                
                // Remove the added constraints
                if (!added_constraints.empty()) {
                    for (Constraint_id cid : added_constraints) {
                        triangulation.cdt.remove_constraint(cid);
                    }
                }
                int polygon_count = triangulation.count_obtuse_triangles();
                if (polygon_count < min_obtuse_triangles) {
                    std::cout<<"Used"<<std::endl;
                    min_obtuse_triangles = polygon_count;
                    improved = true;
                }
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