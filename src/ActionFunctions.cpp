#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include <queue>
#include <cmath>

// Function to check if a triangle is obtuse
bool is_obtuse(const Point& A, const Point& B, const Point& C) {
    
    // Define the angles of the triangle (middle is the vertex)

    // Check if any angle is obtuse
    if (CGAL::angle(A, B, C) == CGAL::OBTUSE || CGAL::angle(B, C, A) == CGAL::OBTUSE || CGAL::angle(C, A, B) == CGAL::OBTUSE) {
        return true;
    }
    return false;
}

// Function to check if two triangles form a convex quadrilateral
bool is_convex_hull(Face_handle fh1, Face_handle fh2) {
    
    // Collect the four points from the two triangles
    std::vector<Point> quadrilateral;

    // Get the vertices of the shared edge
    Point shared_v1 = fh1->vertex(0)->point();
    Point shared_v2 = fh1->vertex(1)->point();
    quadrilateral.push_back(shared_v1);
    quadrilateral.push_back(shared_v2);

    // Add the opposite vertices of the two triangles
    Point tri1_opposite = fh1->vertex(2)->point();
    Point tri2_opposite = fh2->vertex(2)->point();
    quadrilateral.push_back(tri1_opposite);
    quadrilateral.push_back(tri2_opposite);

    // Check if the quadrilateral is convex
    return CGAL::is_convex_2(quadrilateral.begin(), quadrilateral.end());
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

    // Check if the edge is part of the polygon boundary
    Segment_2 edge(v1, v2);
    for (auto it = polygon.edges_begin(); it != polygon.edges_end(); ++it) {
        if (edge == *it) {
            // Edge is part of the boundary, therefore valid
            return true;
        }
    }

    // Check if the edge intersects any polygon edges
    for (auto it = polygon.edges_begin(); it != polygon.edges_end(); ++it) {
        if (CGAL::do_intersect(edge, *it)) {
            // Edge intersects with the polygon boundary
            return false;
        }
    }

    // Edge is valid if no intersections were found and vertices are valid
    return true;
}

// Function to add Steiner points
bool add_optimal_steiner(Triangulation& triangulation) {
    
    // Flag to indicate if the triangulation was improved
    bool improved = false;

    Point best_steiner_point;
    int min_obtuse_triangles = triangulation.count_obtuse_triangles();

    // Store triangluation states
    std::priority_queue<TriangulationState> pq;    
    
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
    
    // For debugging
    std::cout << min_obtuse_triangles << std::endl;

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

        // Calculate squared distances of the edges
        K::FT lAB = CGAL::squared_distance(p0, p1);
        K::FT lBC = CGAL::squared_distance(p1, p2);
        K::FT lCA = CGAL::squared_distance(p2, p0);

        // Calculate max and min edges
        K::FT lmax = std::max({lAB, lBC, lCA});
        K::FT lmin = std::min({lAB, lBC, lCA});

        // Get multiple candidate steiner points positions
        Triangle_2 triangle(obtuse_angle, edge_vertex_1, edge_vertex_2);
        Point circumcenter = CGAL::circumcenter(triangle);
        Point centroid = CGAL::centroid(triangle);
        Point projection = project_point_onto_line(obtuse_angle, edge_vertex_1, edge_vertex_2);
        Point midpoint;
        Point convex_centroid;

        std::vector<Point> candidate_points = {centroid};

        // Calculate the midpoint of the largest edge
        if (lmax == lAB) {
            midpoint = CGAL::midpoint(p0, p1);
        } else if (lmax == lBC) {
            midpoint = CGAL::midpoint(p1, p2);
        } else {
            midpoint = CGAL::midpoint(p2, p0);
        }

        // For each edge of the obtuse triangle check its neighboring triangle
        for (int i = 0; i < 3; ++i) {
            Face_handle neighbor_fh = fh->neighbor(i);

            // Check if the neighbor is finite
            if (!triangulation.cdt.is_infinite(neighbor_fh)) {
                // Check if the obtuse triangle and its neighbor form a convex hull
                if (is_convex_hull(fh, neighbor_fh)) {
                    // Get the vertices of the triangles
                    Point shared_v1 = fh->vertex((i+1) % 3)->point();
                    Point shared_v2 = fh->vertex((i+2) % 3)->point();
                    Point tri1_opposite = fh->vertex(i)->point();
                    Point tri2_opposite = neighbor_fh->vertex(neighbor_fh->index(fh->vertex(i)))->point();

                    // Compute the centroid of the quadrilateral
                    convex_centroid = CGAL::centroid(shared_v1, shared_v2, tri1_opposite, tri2_opposite);
                    
                    candidate_points.push_back(convex_centroid);
                }
            }
        }

        // Validation check to stay inbound
        if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)) {
            candidate_points.push_back(midpoint);
        }
        if (!(triangulation.cdt.is_infinite(triangulation.cdt.locate(circumcenter)))) {
            candidate_points.push_back(circumcenter);
        }
        if (is_edge_valid(edge_vertex_1, edge_vertex_2, triangulation.polygon)) {
            candidate_points.push_back(projection);
        }

        // Iterate through the candidate points
        for (Point& candidate : candidate_points) {
            // Create a copy of the triangulation for testing
            Triangulation tcopy(triangulation);
            tcopy.mark_domain();
            tcopy.cdt.insert(candidate);

            int obtuse_count = tcopy.count_obtuse_triangles();
            
            // For debugging
            std::cout << obtuse_count << std::endl;
            
            // If this approach eliminates obtuse triangles, store it
            if (obtuse_count < min_obtuse_triangles) {
                std::cout << "jere" << std::endl;
                min_obtuse_triangles = obtuse_count;
                best_steiner_point = candidate;
                current_steiner_points.push_back(best_steiner_point);
                pq.push({current_steiner_points, obtuse_count});
                improved = true;
            }
        }

    }

    // Retrieve the state with the least obtuse triangles
    TriangulationState best_state = pq.top();

    // For debugging
    std::cout << "!" << best_state.obtuse_triangle_count << std::endl;

    // If the Steiner point addition improved the triangulation
    if (improved) {
        // Add the Steiner points at the initial triangulation
        for (auto i = 0; i < best_state.steiner_points.size(); ++i) {
            triangulation.cdt.insert(best_state.steiner_points[i]);
            triangulation.polygon.push_back(best_state.steiner_points[i]);
            
            // For debugging
            std::cout << "?" << best_state.steiner_points[i] << std::endl;
        }
    }

    return improved;
}

// Function to eliminate obtuse triangles
void eliminate_obtuse_triangles(Triangulation& triangulation) {
    bool improved = true;

    // Loop terminates when adding steiner points doesn't improve the triangulations
    while (improved) {
        improved = add_optimal_steiner(triangulation);
    }
}