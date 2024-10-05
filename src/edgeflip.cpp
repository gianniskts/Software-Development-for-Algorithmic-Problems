#include "../includes/edgeflip.h"

#include <CGAL/intersections.h>
#include <set>

// Function to calculate angle at a vertex
double angle(const Point& p1, const Point& p2, const Point& p3) {
    // Vectors defined by points
    K::Vector_2 v1 = p2 - p1;
    K::Vector_2 v2 = p3 - p1;

    // Dot product and norms of vectors
    double dot = v1 * v2;
    double norm1_squared = v1 * v1;
    double norm2_squared = v2 * v2;

    // Using Cosine rule return the angle in degrees
    return std::acos(dot / (norm1_squared * norm2_squared)) * 180 / M_PI;
}

// Function to check if we must flip edge
bool should_flip(const CDT& cdt, Edge edge) {
    if (cdt.is_valid()) {
        return false;
    }
    
    // One of the faces sharing the edge
    Face_handle f1 = edge.first;
    // Its opposite face
    Face_handle f2 = f1->neighbor(edge.second);

    // Get the vertices of the triangles involved in the edge flip
    Point A = f1->vertex((edge.second + 1) % 3)->point();
    Point B = f1->vertex((edge.second + 2) % 3)->point();
    Point C = f1->vertex(edge.second)->point();
    Point D = f2->vertex(cdt.mirror_index(edge.first, edge.second))->point(); // Get vertex opposite to the edge in face 2

    double a1 = angle(A, B, C);
    double a2 = angle(A, B, D);

    // Check if the angles are obtuse
    return a1 > 90 || a2 > 90;
}


typedef CDT::Vertex_handle Vertex_handle;
// For debugging purposes
void print_vertex(Vertex_handle v) {
    std::cout << "(" << v->point().x() << ", " << v->point().y() << ")";
}

// Function to perform edge flips
void edge_flip(CDT& cdt) {
    bool flipped;

    //Debugging
    std::set<std::pair<Face_handle, int>> flipped_edges;

    do {
        flipped = false;
        for (auto edge_it = cdt.finite_edges_begin(); edge_it != cdt.finite_edges_end(); ++edge_it) {
            
            // Get the face and edge index from the edge iterator
            Face_handle face = edge_it->first;  // The face containing the edge
            int edge_index = edge_it->second; // The index of the edge in the face

            // Check if the edge is on the boundary
            Face_handle neighbor_face = face->neighbor(edge_index);

            if (cdt.is_infinite(neighbor_face) || cdt.is_infinite(face)) { // Debugging
                std::cout << "Skipping edge on boundary: " 
                          << face->vertex(0)->point() << " "
                          << face->vertex(1)->point() << " "
                          << face->vertex(2)->point() << std::endl;
                
                continue; // Skip this edge if it's on the boundary
            }

            // Debugging
            if (face == neighbor_face) {
                std::cout << "Skipping flip for edge with same face: ";
                print_vertex(face->vertex(0)); std::cout << " ";
                print_vertex(face->vertex(1)); std::cout << " ";
                print_vertex(face->vertex(2)); std::cout << std::endl;
                continue; // Avoid flipping if the two faces are the same
            }

            // Perform the flip
            if (should_flip(cdt, *edge_it)) {
                // Check if this edge has already been flipped
                if (flipped_edges.find(std::make_pair(face, edge_index)) != flipped_edges.end()) {
                    std::cout << "Skipping already flipped edge: ";
                    print_vertex(face->vertex(0)); std::cout << " ";
                    print_vertex(face->vertex(1)); std::cout << " ";
                    print_vertex(face->vertex(2)); std::cout << std::endl;
                    continue;
                }

                std::cout << "Flipping edge: ";
                print_vertex(face->vertex(0)); std::cout << " ";
                print_vertex(face->vertex(1)); std::cout << " ";
                print_vertex(face->vertex(2)); std::cout << std::endl;

                // Flip the edge
                cdt.tds().flip(face, edge_index);
                flipped_edges.insert(std::make_pair(face, edge_index)); // Mark as flipped
                flipped = true;
            }
        }
    } while (flipped); // Repeat until no more flips are possible
}