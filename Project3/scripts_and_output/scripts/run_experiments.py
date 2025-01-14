import csv
import glob
import itertools
import subprocess
import json
import time
import os

# IMPORTANT: MOVE THE SCRIPT TO THE ROOT DIRECTORY OF THE PROJECT

def run_one_experiment(instance_path, output_path, extra_args):
    """
    Runs one instance of ./build/main with the given arguments.
    Returns (exitcode, elapsed_time_seconds).
    """
    cmd = [
        "./build/main",
        "-i", instance_path,
        "-o", output_path
    ] + extra_args

    t0 = time.time()
    try:
        process = subprocess.run(cmd, capture_output=True, timeout=20)  # Set a timeout of 20 seconds
        t1 = time.time()
        elapsed = t1 - t0

        exitcode = process.returncode
        if exitcode != 0:
            print("[ERROR] cmd:", cmd)
            print("stderr:", process.stderr.decode("utf-8"))
            print("stdout:", process.stdout.decode("utf-8"))

        return exitcode, elapsed
    except subprocess.TimeoutExpired:
        print(f"[TIMEOUT] Command exceeded 20 seconds: {cmd}")
        return -1, 20  # Return -1 as exit code and 20 seconds as elapsed time


def parse_output_json(output_path):
    """
    Reads the output JSON from 'output_path' and extracts relevant metrics:
    - method
    - obtuse_count
    - number of Steiner points
    - randomization
    - etc.
    Returns a dictionary with these values.
    """
    if not os.path.exists(output_path):
        return {}

    with open(output_path, "r") as f:
        data = json.load(f)
        result = {
            "method": data.get("method", ""),
            "obtuse_count": data.get("obtuse_count", -1),
            "randomization": data.get("randomization", False),
            "energy": data.get("energy", -99.9),
            "p_bar": data.get("p_bar", -99.9),
            "category": data.get("category", "default")
        }
        # For number_of_steiner_points:
        spx = data.get("steiner_points_x", [])
        result["steiner_points"] = len(spx)

        return result


def main():
    instance_folder = "challenge_instances_cgshop25"
    os.makedirs("best_outputs", exist_ok=True)

    instances = glob.glob(os.path.join(instance_folder, "*.json"))

    if not instances:
        print(f"[ERROR] No .json files found in {instance_folder}")
        return
    
    start_test = 1
    end_test = 150
    instances = instances[start_test - 1:end_test]
    instances_len = len(instances)

    methods = ["local", "sa", "ant"]
    steiner_methods = ["projection", "midpoint", "centroid", "circumcenter", "mean_of_adjacent"]

    # Generate all combinations of Steiner methods
    steiner_combinations = []
    for i in range(1, len(steiner_methods) + 1):
        steiner_combinations.extend(itertools.combinations(steiner_methods, i))

    csv_file = "final_results_FINAL.csv"
    best_results_file = "best_results_FINAL.csv"
    best_output_file = "best_output_FINAL.csv"

    file_exists = os.path.isfile(csv_file)
    best_file_exists = os.path.isfile(best_results_file)
    best_output_exists = os.path.isfile(best_output_file)

    instance_counter = 0
    total_tests = len(instances) * len(methods) * len(steiner_combinations) * 2

    with open(csv_file, "a", newline="") as f, \
         open(best_results_file, "a", newline="") as best_f, \
         open(best_output_file, "a", newline="") as best_out_f:

        writer = csv.DictWriter(f, fieldnames=[
            "instance", "method", "steiner_methods",
            "randomization", "status", "time",
            "obtuse_count", "steiner_points", "energy", "p_bar", "category"
        ])
        best_writer = csv.DictWriter(best_f, fieldnames=[
            "instance", "method", "steiner_methods",
            "randomization", "status", "time",
            "obtuse_count", "steiner_points", "energy", "p_bar", "category"
        ])
        best_output_writer = csv.DictWriter(best_out_f, fieldnames=[
            "instance", "output_json"
        ])

        if not file_exists:
            writer.writeheader()
        if not best_file_exists:
            best_writer.writeheader()
        if not best_output_exists:
            best_output_writer.writeheader()

        for instance_path in instances:
            best_result = None
            best_result_output_file = None
            instance_counter += 1
            print(f"Running instance {instance_counter}/{instances_len}: {instance_path}")
            basename = os.path.basename(instance_path)

            # Track consecutive failures
            consecutive_failures = 0

            for method_name in methods:
                # If we've already hit 5 consecutive fails, break to next instance
                if consecutive_failures >= 5:
                    print(f"[INFO] 5 consecutive failures. Skipping remaining tests for {basename}.")
                    break

                for steiner_set in steiner_combinations:
                    # If we've already hit 5 consecutive fails, break
                    if consecutive_failures >= 5:
                        print(f"[INFO] 5 consecutive failures. Skipping remaining tests for {basename}.")
                        break

                    chosen_steiner_methods = ",".join(steiner_set)

                    for random_flag in [True, False]:
                        # If we've already hit 5 consecutive fails, break
                        if consecutive_failures >= 5:
                            print(f"[INFO] 5 consecutive failures. Skipping remaining tests for {basename}.")
                            break

                        rand_suffix = "rand" if random_flag else "noRand"
                        output_file = f"results/{basename}__{method_name}__{rand_suffix}__{chosen_steiner_methods}.json"

                        extra_args = [
                            "-methods", method_name,
                            "-steiner", chosen_steiner_methods
                        ]
                        if not random_flag:
                            extra_args.append("--no-randomization")

                        print(
                            f"Currently testing: method={method_name}, "
                            f"steiner={chosen_steiner_methods}, "
                            f"randomization={random_flag}, "
                            f"instance={instance_counter}/{instances_len}"
                        )

                        exitcode, elapsed_sec = run_one_experiment(
                            instance_path, output_file, extra_args
                        )

                        if exitcode != 0:
                            # Failed experiment, increment consecutive fail count
                            consecutive_failures += 1
                            writer.writerow({
                                "instance": basename,
                                "method": method_name,
                                "steiner_methods": chosen_steiner_methods,
                                "randomization": random_flag,
                                "status": "FAILED",
                                "time": elapsed_sec,
                                "obtuse_count": -1,
                                "steiner_points": -1,
                                "energy": -99.9,
                                "p_bar": -99.9,
                                "category": "default"
                            })

                            # If we reached 5 consecutive fails, break out
                            if consecutive_failures >= 5:
                                print(f"[INFO] Reached 5 consecutive failures for {basename}. Moving to next instance.")
                                break
                            continue
                        else:
                            # Success, reset consecutive failures
                            consecutive_failures = 0

                        # Parse the JSON and write results
                        data = parse_output_json(output_file)
                        writer.writerow({
                            "instance": basename,
                            "method": data.get("method", method_name),
                            "steiner_methods": chosen_steiner_methods,
                            "randomization": data.get("randomization", random_flag),
                            "status": "OK",
                            "time": elapsed_sec,
                            "obtuse_count": data.get("obtuse_count", -1),
                            "steiner_points": data.get("steiner_points", -1),
                            "energy": data.get("energy", -99.9),
                            "p_bar": data.get("p_bar", -99.9),
                            "category": data.get("category", "default")
                        })

                        # Update best result
                        p_bar = data.get("p_bar", -99.9)
                        steiner_points = data.get("steiner_points", float('inf'))
                        energy = data.get("energy", float('inf'))

                        if p_bar > 0:  # Converged
                            if (best_result is None or
                                    p_bar > best_result["p_bar"] or
                                    (p_bar == best_result["p_bar"] and steiner_points < best_result["steiner_points"])):
                                best_result = {
                                    **data,
                                    "instance": basename,
                                    "status": "BEST",
                                    "time": elapsed_sec
                                }
                                best_result_output_file = output_file
                        else:  # Did not converge
                            # If best_result is None or we only have non-convergent results so far,
                            # pick the one with the lower energy (i.e., better solution):
                            if (best_result is None or
                                    (best_result["p_bar"] <= 0 and energy < best_result["energy"])):
                                best_result = {
                                    **data,
                                    "instance": basename,
                                    "status": "BEST",
                                    "time": elapsed_sec
                                }
                                best_result_output_file = output_file

            if best_result:
                # Write to best_results.csv
                best_writer.writerow(best_result)

                # Also store the entire output JSON in best_output.csv
                if best_result_output_file and os.path.isfile(best_result_output_file):
                    with open(best_result_output_file, "r") as json_f:
                        entire_json_str = json_f.read()
                        best_output_writer.writerow({
                            "instance": basename,
                            "output_json": entire_json_str
                        })

    print("Done! Results saved in final_results_NEW.csv, best_results.csv, and best_output.csv.")


if __name__ == "__main__":
    main()
