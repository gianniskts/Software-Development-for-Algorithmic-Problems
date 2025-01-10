#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include "../includes/Parsing.h"
#include <CGAL/convex_hull_2.h>
#include <queue>
#include <cmath>
#include <random>
#include <chrono>

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

/* Roll-back to 1st Project implementation (if delaunay boolean parameter is true)*/

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
        Point convex_centroid;

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
void eliminate_obtuse_triangles(Triangulation& triangulation) {
    bool improved = true;

    // Loop terminates when adding steiner points doesn't improve the triangulations
    while (improved) {
        improved = add_optimal_steiner(triangulation);
    }

}

// Function to generate all subsets of a vector of Steiner sub-methods
// Each subset is a vector<string>
std::vector<std::vector<std::string>> generate_subsets(const std::vector<std::string> &items)
{
    std::vector<std::vector<std::string>> subsets;
    size_t n = items.size();
    // We'll generate 2^n subsets
    for (size_t mask = 1; mask < (1UL << n); ++mask)
    {
        std::vector<std::string> subset;
        for (size_t bit = 0; bit < n; ++bit)
        {
            if (mask & (1UL << bit))
            {
                subset.push_back(items[bit]);
            }
        }
        subsets.push_back(subset);
    }
    return subsets;
}

static bool is_polygon_convex(const Polygon_2 &poly)
{
    return poly.is_convex();
}

// Check if all edges of the given boundary polygon are strictly horizontal or vertical.
static bool all_boundary_edges_axis_aligned(const Polygon_2 &boundary)
{
    for (std::size_t i = 0; i < boundary.size(); i++)
    {
        const Point &p1 = boundary[i];
        const Point &p2 = boundary[(i + 1) % boundary.size()];
        double x1 = CGAL::to_double(p1.x());
        double y1 = CGAL::to_double(p1.y());
        double x2 = CGAL::to_double(p2.x());
        double y2 = CGAL::to_double(p2.y());
        // If neither x nor y is the same, it's not axis-aligned
        if (x1 != x2 && y1 != y2)
        {
            return false;
        }
    }
    return true;
}

// Function to check if the boundary is identical to the convex envelope
static bool is_identical_to_convex_envelope(const std::vector<Point> &boundary,
                                            const std::vector<Point> &all_points)
{
    // Step 1: Reconstruct the convex hull of all points
    std::vector<Point> convex_hull_points;
    CGAL::convex_hull_2(all_points.begin(), all_points.end(), std::back_inserter(convex_hull_points));

    // Step 2: Compare the convex hull with the boundary
    // Ensure both have the same size
    if (convex_hull_points.size() != boundary.size())
    {
        return false;
    }

    // Step 3: Ensure the boundary points match the convex hull points in the same order
    // Normalize both point sequences to start from the same point for comparison
    auto normalize_sequence = [](std::vector<Point> &points)
    {
        auto min_it = std::min_element(points.begin(), points.end());
        std::rotate(points.begin(), min_it, points.end());
    };

    std::vector<Point> normalized_boundary = boundary;
    normalize_sequence(convex_hull_points);
    normalize_sequence(normalized_boundary);

    // Step 4: Compare normalized sequences
    return std::equal(convex_hull_points.begin(), convex_hull_points.end(), normalized_boundary.begin());
}

// Check to detect "closed" constraints: any constraint lying exactly on the boundary
// or forming a cycle in the interior. If true, we classify as "closed".
static bool is_constraints_closed_or_on_boundary(const InputJSON &input,
                                                 const Polygon_2 &boundary)
{
    // 1) Make a set of boundary edges in canonical form (minIndex,maxIndex)
    std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> boundaryEdges;
    auto makeEdge = [&](int a, int b)
    {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };
    for (std::size_t i = 0; i < input.region_boundary.size(); i++)
    {
        int curr = input.region_boundary[i];
        int nxt = input.region_boundary[(i + 1) % input.region_boundary.size()];
        boundaryEdges.insert(makeEdge(curr, nxt));
    }

    // 2) Check if any constraint is exactly a boundary edge
    for (auto &ac : input.additional_constraints)
    {
        auto e = makeEdge(ac.first, ac.second);
        if (boundaryEdges.find(e) != boundaryEdges.end())
        {
            // This constraint is exactly on the boundary --> "closed"
            return true;
        }
    }

    // 3) Build adjacency to detect interior cycles
    std::unordered_map<int, std::vector<int>> adj;
    for (auto &ac : input.additional_constraints)
    {
        adj[ac.first].push_back(ac.second);
        adj[ac.second].push_back(ac.first);
    }

    std::unordered_set<int> visited;
    std::function<bool(int, int)> dfsDetectCycle = [&](int node, int parent) -> bool
    {
        visited.insert(node);
        for (int nbr : adj[node])
        {
            if (nbr == parent)
                continue; // skip edge we came from
            if (visited.find(nbr) != visited.end())
            {
                // Found a back edge => cycle
                return true;
            }
            else if (dfsDetectCycle(nbr, node))
            {
                return true;
            }
        }
        return false;
    };

    // 4) DFS from each unvisited node to check if there's a cycle
    for (auto &kv : adj)
    {
        int start = kv.first;
        if (visited.find(start) == visited.end())
        {
            if (dfsDetectCycle(start, -1))
            {
                return true; // found a cycle => "closed"
            }
        }
    }

    // Otherwise, constraints are "open"
    return false;
}

std::string detect_category(const InputJSON &input)
{
    // 1) Reconstruct the region boundary polygon
    Polygon_2 boundary;
    for (int idx : input.region_boundary)
    {
        boundary.push_back(Point(input.points_x[idx], input.points_y[idx]));
    }

    // Convert boundary polygon to a vector of Points
    std::vector<Point> boundary_points;
    for (auto it = boundary.vertices_begin(); it != boundary.vertices_end(); ++it)
    {
        boundary_points.push_back(*it);
    }

    // Convert all input points to a vector of CGAL Points for convex envelope verification
    std::vector<Point> all_points;
    for (size_t i = 0; i < input.points_x.size(); ++i)
    {
        all_points.push_back(Point(input.points_x[i], input.points_y[i]));
    }

    // 2) Check if boundary is convex
    bool polygon_is_convex = is_polygon_convex(boundary);
    int c = input.num_constraints;

    // 3) Handle the convex vs. non-convex boundary separately
    if (polygon_is_convex)
    {
        bool is_convex_envelope = is_identical_to_convex_envelope(boundary_points, all_points);
        if (!is_convex_envelope)
        {
            return "Not Convex Envelope";
        }
        // --- Category A ---
        // Convex boundary, no constraints
        if (c == 0)
        {
            return "A";
        }

        // We have constraints -> decide if "closed" or "open"
        bool closedConstraints = is_constraints_closed_or_on_boundary(input, boundary);

        if (closedConstraints)
        {
            // --- Category C ---
            // Convex boundary with constraints that form closed polygons or lie on boundary
            return "C";
        }
        else
        {
            // --- Category B ---
            // Convex boundary, constraints are "open"
            return "B";
        }
    }
    else
    {
        // Non-convex boundary -> check for Category D or else E
        // Category D: Non-convex boundary with axis-aligned edges and no constraints
        bool boundary_is_axis_aligned = all_boundary_edges_axis_aligned(boundary);
        if (c == 0 && boundary_is_axis_aligned)
        {
            return "D";
        }
        // Otherwise -> E
        return "E";
    }
}


// Inserting one random Steiner point near the centroid of a randomly chosen
// obtuse triangle, using a Gaussian distribution for the (dx, dy) offset
// from the exact centroid. Param sigma is the standard deviation for the 2D normal distribution
bool insert_gaussian_near_centroid_of_obtuse_triangle(
    Triangulation &triang,
    double sigma
) {
    // 1) Collect all obtuse triangles in the domain
    std::vector<Face_handle> obtuse_faces;
    for (auto face_it = triang.cdt.finite_faces_begin();
         face_it != triang.cdt.finite_faces_end(); ++face_it)
    {
        if (!triang.is_face_in_domain(face_it)) continue;
        Point p0 = face_it->vertex(0)->point();
        Point p1 = face_it->vertex(1)->point();
        Point p2 = face_it->vertex(2)->point();
        if (is_obtuse(p0, p1, p2)) {
            obtuse_faces.push_back(face_it);
        }
    }

    // If no obtuse triangles, nothing to do
    if (obtuse_faces.empty()) {
        return false;
    }

    // 2) Randomly pick one obtuse face
    static std::mt19937 rng(
        std::chrono::steady_clock::now().time_since_epoch().count()
    );
    std::uniform_int_distribution<size_t> face_dist(0, obtuse_faces.size() - 1);
    Face_handle fh = obtuse_faces[face_dist(rng)];

    // 3) Compute the centroid of this triangle
    Point p0 = fh->vertex(0)->point();
    Point p1 = fh->vertex(1)->point();
    Point p2 = fh->vertex(2)->point();
    Point centroid = CGAL::centroid(p0, p1, p2);

    double cx = CGAL::to_double(centroid.x());
    double cy = CGAL::to_double(centroid.y());

    // 4) Prepare to generate random offsets using a Gaussian distribution
    std::normal_distribution<double> normal_dist(0.0, sigma);

    // 5) Attempt multiple times to get a valid point that is inside the domain
    const int maxAttempts = 10;  // Adjust as needed
    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        double dx = normal_dist(rng);
        double dy = normal_dist(rng);

        Point new_steiner(cx + dx, cy + dy);

        // Check if this point is located inside the domain
        CDT::Locate_type lt;
        int li;
        Face_handle locatedFace = triang.cdt.locate(new_steiner, lt, li);

        // We consider a point valid if:
        //   1) it is not located at an infinite face/edge/vertex
        //   2) the face we found is inside the polygon domain
        if (lt == CDT::FACE && !triang.cdt.is_infinite(locatedFace) && 
            triang.is_face_in_domain(locatedFace))
        {
            // Insert into the CDT
            triang.cdt.insert(new_steiner);
            triang.mark_domain();
            triang.randomizationUsed = true;
            return true;
        }
    }

    // If we exhaust attempts without finding an in-domain point, return false
    return false;
}


// p(n) computation:  |log(obtuse_prev/obtuse_curr)| / |log(n_prev/n_curr)|
double compute_pn(int obtuse_prev, int obtuse_curr, int n_prev, int n_curr) {
    // We add protective checks. Also using absolute value as required.
    if (obtuse_prev <= 0 || obtuse_curr <= 0) {
        // Once we’re at 0 obtuse, we typically stop measuring or return 0
        return 0.0;
    }
    if (n_curr <= 1 || n_prev <= 1) {
        return 0.0;
    }

    double top    = std::log( (double)obtuse_prev / (double)obtuse_curr );
    double bottom = std::log( (double)n_prev     / (double)n_curr     );
    if (std::fabs(bottom) < 1e-12) {
        return 0.0;
    }
    return std::fabs(top) / std::fabs(bottom);
}

void randomize_if_stuck(
    Triangulation &triang, 
    int stepsWithoutImprovement,
    int threshold,
    bool randomization_enabled
) {
    if (!randomization_enabled) {
        return;
    }

    static int randomCountGlobal = 0;
    static const int maxRandomGlobal = 3;

    if (stepsWithoutImprovement >= threshold && (randomCountGlobal < maxRandomGlobal))
    {
        std::cout << "[Randomization Triggered] Inserting random Steiner points near centroids...\n";

        // Example: We'll do 3 random inserts
        // Each insert is near the centroid of a random obtuse triangle
        for (int i = 0; i < 3; i++) {
            bool success = insert_gaussian_near_centroid_of_obtuse_triangle(triang, /* sigma= */ 50.0);
            // If we found no obtuse triangle or insertion failed, we can break
            if (!success) break;
        }

        ++randomCountGlobal;
    }
}

