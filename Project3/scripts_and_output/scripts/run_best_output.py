import os
import glob
import subprocess

# IMPORTANT: MOVE THE SCRIPT TO THE ROOT DIRECTORY OF THE PROJECT

def main():
    instance_folder = "challenge_instances_cgshop25"
    output_folder = "best_outputs2"
    
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Get all .json files in the instance folder
    instances = glob.glob(os.path.join(instance_folder, "*.json"))
    print(f"Found {len(instances)} .json files in {instance_folder}")

    if not instances:
        print(f"[ERROR] No .json files found in {instance_folder}")
        return

    # Iterate through each instance and process it
    for instance_path in instances:
        # Generate output path
        instance_name = os.path.basename(instance_path).replace(".json", ".output.json")
        output_path = os.path.join(output_folder, instance_name)

        # Command to run the C++ program
        cmd = [
            "./build/main",
            "-i", instance_path,
            "-o", output_path,
            "-preselected_params"
        ]

        try:
            # Execute the command with a timeout of 60 seconds
            print(f"Processing: {instance_path}")
            subprocess.run(cmd, timeout=180, check=True)
            print(f"Successfully processed {instance_path}, output saved to {output_path}")
        except subprocess.TimeoutExpired:
            print(f"[ERROR] Timeout occurred while processing {instance_path}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] An error occurred while processing {instance_path}: {e}")

if __name__ == "__main__":
    main()
