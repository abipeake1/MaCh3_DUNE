import yaml
import os

def update_pot_values(directory, total_pot, custom_pot_file, custom_pot_percentage):
    yaml_files = [f for f in os.listdir(directory) if f.endswith(".yaml")]
    
    if custom_pot_file not in yaml_files:
        print(f"Error: {custom_pot_file} not found in {directory}")
        return
    
    # Number of files excluding the one getting a custom POT
    other_files = [f for f in yaml_files if f != custom_pot_file]
    num_other_files = len(other_files)
    
    if num_other_files == 0:
        print("Error: No other YAML files found to distribute POT.")
        return
    
    custom_pot_value = total_pot * (custom_pot_percentage / 100)
    remaining_pot = total_pot - custom_pot_value
    pot_per_file = remaining_pot / num_other_files
    
    # Function to update a YAML file
    def modify_yaml(file_path, new_pot):
        with open(file_path, 'r') as f:
            data = yaml.safe_load(f)
        
        if "POT" in data:
            data["POT"] = new_pot
        else:
            print(f"Warning: No POT key found in {file_path}, skipping.")
            return
        
        with open(file_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False)
    
    # Update the custom POT file
    modify_yaml(os.path.join(directory, custom_pot_file), custom_pot_value)
    
    # Update all other files
    for file in other_files:
        modify_yaml(os.path.join(directory, file), pot_per_file)
    
    print("POT values updated successfully.")

# Example usage
update_pot_values("./configs/Samples/OA_Samples_subsamples", 7.238e21, "0m_all.yaml", 50)
