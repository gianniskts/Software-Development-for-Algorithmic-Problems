#ifndef POLYGONMANIPULATION_H
#define POLYGONMANIPULATION_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Segment_2.h>
#include <iostream>
#include <unordered_map>
#include <boost/property_map/property_map.hpp>

template <typename T> struct InputJSON;

// Necessary declarations
typedef CGAL::Exact_predicates_exact_constructions_kernel           K;
typedef CGAL::Triangulation_vertex_base_2<K>                        Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>              Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                 TDS;
typedef CGAL::Exact_predicates_tag                                  Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>    CDT;
typedef CDT::Face_handle                                            Face_handle;
typedef CDT::Point                                                  Point;
typedef CDT::Edge                                                   Edge;
typedef CDT::Face_iterator                                          Face_iterator;
typedef CGAL::Polygon_2<K>                                          Polygon_2;
typedef CGAL::Triangle_2<K>                                         Triangle_2;
typedef CGAL::Segment_2<K>                                          Segment_2;

// Function to perform constrained delaunay triangulation
template <typename T>
Polygon_2 delaunay_const_triangulation(InputJSON<T> input_data);

#endif // POLYGONMANIPULATION_H