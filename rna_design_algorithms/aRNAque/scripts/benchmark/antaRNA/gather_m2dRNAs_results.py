import subprocess
import pandas as pd

df = pd.read_csv("Eterna100.txt")


structures = df["str"].to_list()

# Dictionary to store the results
results = {}

# Loop through each structure and run the command
counter = 0
for structure in structures:
    counter += 1
    # command = f'./m2dRNAs 1 "{structure}" 52 50 1.0'
    command = f'./m2dRNAs 1 "{structure}" 50 -1 0 > $output-{counter}.out'
    try:
        # Execute the command and capture the output
        result = subprocess.check_output(command, shell=True, universal_newlines=True)

        # Store the output in the dictionary
        results[structure] = result
    except subprocess.CalledProcessError as e:
        # Handle errors in command execution
        print(f"Error executing command for structure {structure}: {e}")
        results[structure] = f"Error: {e}"