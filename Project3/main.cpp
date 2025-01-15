#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <sstream>
#include "../includes/Parsing.h"
#include "../includes/PolygonManipulation.h"
#include "../includes/ActionFunctions.h"
#include "../includes/MyTriangulation.h"
#include "../includes/OptimizationMethods.h"
#include "../includes/draw.h"

int main(int argc, char *argv[])
{
    try
    {
        if (argc < 5)
        {
            std::cerr << "Usage: ./opt_triangulation -i /path/to/input.json "
                         "-o /path/to/output.json "
                         "[options...]\n\n"
                      << "Options:\n"
                      << "  -preselected_params           Use the 'auto' code path with preselected param combos.\n"
                      << "  -methods local,sa,ant         Comma-separated list of methods to run (only local|sa|ant).\n"
                      << "  -steiner projection,mid,...   Comma-separated list of Steiner sub-methods to use (projection|mid|centroid|circumcenter|adjacent). If omitted, use all.\n"
                      << "  -test_subsets                 If provided, test all subsets of the chosen Steiner sub-methods.\n"
                      << "  --no-randomization            Disable random insertion of Steiner points during deadlocks.\n"
                      << "  -advance_output               Adds p_bar, energy in the output file"
                      << std::endl;
            return 1;
        }

        std::string input_file;
        std::string output_file;

        // Flags
        bool use_preselected = false;         // Use the preselected method with the best parameters
        bool test_all_subsets = false;       // Test all subsets of the chosen Steiner sub-methods
        bool randomization_enabled = true;  // Enable random insertion of Steiner points during deadlocks
        bool advancedOutput = false;       // Adds p_bar, energy in the output file

        // Which methods to run
        std::vector<std::string> chosen_methods;
        // Which Steiner techniques to use
        std::vector<std::string> chosen_steiner_methods;

        // Command-line argument handling
        for (int i = 1; i < argc; i++)
        {
            std::string arg = argv[i];
            if ((arg == "-i" || arg == "--input") && i + 1 < argc)
            {
                input_file = argv[++i];
            }
            else if ((arg == "-o" || arg == "--output") && i + 1 < argc)
            {
                output_file = argv[++i];
            }
            else if (arg == "-preselected_params")
            {
                use_preselected = true;
            }
            else if (arg == "-methods" && i + 1 < argc)
            {
                chosen_methods = split_csv(argv[++i]);
            }
            else if (arg == "-steiner" && i + 1 < argc)
            {
                chosen_steiner_methods = split_csv(argv[++i]);
            }
            else if (arg == "-test_subsets")
            {
                test_all_subsets = true;
            }
            else if (arg == "--no-randomization")
            {
                randomization_enabled = false;
            }
            else if (arg == "-advanced_output")
            {
                advancedOutput = true;
            }
            else
            {
                std::cerr << "Unknown or incomplete argument: " << arg << std::endl;
                return 1;
            }
        }

        // Parse the input file
        InputJSON input_json = parse_file(input_file);
        input_json.randomization_enabled = randomization_enabled;
        if (input_json.method == "auto") 
            use_preselected = true;
        
        if (use_preselected) // Use the preselected method with the best parameters
        {
            chosen_methods.push_back("sa");
            std::string category = detect_category(input_json);
            // std::cout << "[auto_method] Detected category: " << category << std::endl;
            if (category == "A")
            {
                // A => SA + {"projection"}
                chosen_steiner_methods = {"projection"};
            }
            else if (category == "B")
            {
                // B => SA + {"projection", "circumcenter"}
                chosen_steiner_methods = {"projection", "circumcenter"};
            }
            else if (category == "C")
            {
                // C => SA + {"projection", "circumcenter", "centroid"}
                chosen_steiner_methods = {"projection", "circumcenter", "centroid"};
            }
            else if (category == "D")
            {
                // D => SA + {"projection", "circumcenter", "midpoint"}
                chosen_steiner_methods = {"projection", "circumcenter", "midpoint"};
            }
            else if (category == "E")
            {
                // E => SA + {"projection", "circumcenter"}
                chosen_steiner_methods = {"projection", "circumcenter"};
            }
        }

        // If we have NOT specified methods from the command line,
        // fallback to the method from the JSON (which must be "local", "sa", or "ant").
        // Input with parameters is also supported.
        if (chosen_methods.empty())
        {
            // If user never specified -methods, fallback to JSON
            if (input_json.method == "local" ||
                input_json.method == "sa" ||
                input_json.method == "ant")
            {
                chosen_methods.push_back(input_json.method);
            }
            else
            {
                // By default let's do local if the JSON is missing or invalid
                chosen_methods.push_back("local");
            }
        }

        // If no Steiner subset specified, we default to "projection,mid,centroid,circumcenter,mean_of_adjacent"
        if (chosen_steiner_methods.empty())
        {
            std::cout << "[INFO] No Steiner methods specified. Using default set.\n";
            chosen_steiner_methods = {"projection", "midpoint", "centroid", "circumcenter", "mean_of_adjacent"};
        }

        // We'll decide which subsets to use
        std::vector<std::vector<std::string>> all_steiner_subsets;
        if (test_all_subsets)
        {
            all_steiner_subsets = generate_subsets(chosen_steiner_methods);
        }
        else
        {
            // Just one "subset", which is the entire chosen set
            all_steiner_subsets.push_back(chosen_steiner_methods);
            // std::cout << "Using all Steiner methods: " << vector_to_string(chosen_steiner_methods) << std::endl;
        }

        // Helper to compute the "energy" of a triangulation for tie-breaking or picking best.
        // IMPORTANT: pass the Triangulation *by value* so we're not capturing references
        // to out-of-scope data inside the lambda.
        auto compute_energy = [&](Triangulation t)
        {
            // Safety check: if no faces, we skip or return a large penalty (so it won't be best).
            if (!t.cdt.is_valid() || t.cdt.number_of_faces() == 0)
            {
                std::cerr << "[WARNING] Triangulation is invalid or empty. "
                             "Cannot count obtuse triangles.\n";
                return 1e15;
            }

            std::cout << "Computing energy for triangulation...\n";
            // Now safe to call count_obtuse_triangles
            int obtuse = t.count_obtuse_triangles();
            t.mark_domain();
            // IF WE WANT TO DRAW THE TRIANGULATION
            CGAL::draw(t.cdt, t.in_domain);

            // Steiner points are the difference between the number of vertices and the number of input points
            int steiner = t.cdt.number_of_vertices() - input_json.num_points;
            return input_json.alpha * obtuse + input_json.beta * steiner;
        };

        Triangulation best_overall;           // store the best triangulation
        double best_energy = 1e15;           // store the best energy
        bool have_best = false;             // flag to check if we have a valid triangulation
        InputJSON best_input = input_json; // store the best input

        // For each chosen method, run it on each subset of Steiner methods
        for (const auto &method_name : chosen_methods)
        {
            // std::cout << "Running method: " << method_name << std::endl;
            for (const auto &steiner_subset : all_steiner_subsets)
            {
                // std::cout << "  Subset: " << vector_to_string(steiner_subset) << std::endl;

                // Copy input JSON but override with the chosen Steiner subset
                InputJSON local_input = input_json;
                local_input.method = method_name;
                local_input.steiner_methods = steiner_subset;

                // Produce a triangulation
                Triangulation T;
                if (method_name == "local")
                {
                    T = local_search(local_input);
                }
                else if (method_name == "sa")
                {
                    T = simulated_annealing(local_input);
                }
                else if (method_name == "ant")
                {
                    std::cout << "Running ant_colony_optimization with "
                              << local_input.kappa << " ants.\n";
                    T = ant_colony_optimization(local_input);
                }
                else
                {
                    throw std::runtime_error("Unknown method specified (should never happen).");
                }

                // Ensure T is valid before computing energy
                if (!T.cdt.is_valid())
                {
                    std::cerr << "[WARNING] T returned from " << method_name
                              << " is invalid. Skipping.\n";
                    continue;
                }

                // Now compute the energy using T
                double e = compute_energy(T);
                if (e < best_energy)
                {
                    best_energy = e;  // store the best energy
                    best_overall = T; // store the best triangulation
                    have_best = true; // flag to check if we have a valid triangulation
                    best_input = local_input; // store the best input
                }
            }
        }

        // If no valid triangulation was found, do a fallback
        if (!have_best) 
        {
            std::cerr << "[WARNING] No triangulation produced. Returning a Delaunay fallback.\n";
            best_overall = delaunay_const_triangulation(input_json);
        }

        output_results(output_file, best_input, best_overall, advancedOutput);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
