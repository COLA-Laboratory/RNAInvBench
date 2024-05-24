import subprocess
import pandas as pd

def run_mcts(target_structure=None, input_data_dir=None, gc_content=0.5,
             gc_content_dev=0.2, pseudoknots=0, timeout=600):
    """
    run_mcts
    --------
    Parameters:
        target_structure:   The Secondary Structure that we are aiming to find a
                            sequence for. (Formatted with dot bracket notation)
        input_data_dir:     File pathway of the input data (if we have any)
        gc_content:         The target GC-content of the RNA sequence, choose value
                            from the range [0,1].
        gc_content_dev:     The GC-content deviation of the solution, which is in
                            range [0,0.02]. Generally, the lower this value, the more
                            accurate the solution is.
        pseudoknots:        0 if there are no pseudoknots present, 1 if pseudoknots
                            are present.
        timeout:            The maximum amount of time that the the MCTS algorithm
                            can spend on each sequence.

    Returns:

        Search Length:
                            The maximimum amount of iterations the algorithm takes to
                            explore the search space before returning a solution.
        Sequence:
                            The RNA Sequence found by the algorithm.
        Run Time:
                            How long the algorithm took to find the sequence
        GC content:
                            The ratio of G and C nucleotides to the total amount of
                            nucleotides.
        GC distance:
                            The absolute difference between the GC content produced by
                            the algorithm, and the actual GC content in the sequence.
        Normalised Hamming Distance:
                            The Normalised hamming distance from the predicted
                            structure when compared with the real structure.
    """
    if target_structure == None:
        print("Running MCTS for file pathway of sequences...")
        outputs = []
        # Run MCTS for all sequences
        df = pd.read_csv(input_data_dir)
        for index, row in df.iterrows():
            print("Sequence: ", index)
            args = [
                "-s", str(row["str"]),
                "-GC", str(gc_content),
                "-d", str(gc_content_dev),
                "-pk", str(pseudoknots),
                "-to", str(timeout)
            ]
            # Construct the command to execute learna.py with arguments
            command = ["python", "bin/MCTS-RNA.py"] + args

            try:
                output = subprocess.check_output(command)
                output_string = output.decode().replace("\\n", "\n")

            except subprocess.CalledProcessError as e:
                output_string = e.output.decode().replace("\\n", "\n")
                print("Error occurred while executing:", output_string)
                error = True
            print(output_string)
            outputs.append(output_string)
        results = []
        for lines in outputs:
            temp = {}
            line = lines.split("\n")
            for part in line:
                if part.strip():
                    # Split the line by colon
                    parts = part.split(":")
                    # Extract the key and value, removing leading/trailing whitespace
                    key = parts[0].strip()
                    value = parts[1].strip()
                    # Append key-value pair to the dictionary

                    temp[key] = value
            results.append(temp)
        return results
    else:
        # Run MCTS for just one sequence
        args = [
            "-s", str(target_structure),
            "-GC", str(gc_content),
            "-d", str(gc_content_dev),
            "-pk", str(pseudoknots)
        ]
        # Construct the command to execute learna.py with arguments
        command = ["python", "bin/MCTS-RNA.py"] + args

        print("Running MCTS...")
        try:
            output = subprocess.check_output(command)
            output_string = output.decode().replace("\\n", "\n")

        except subprocess.CalledProcessError as e:
            output_string = e.output.decode().replace("\\n", "\n")
            print("Error occurred while executing:", output_string)
            error = True
        results = []
        lines = output_string.split("\n")
        print("lines: ", lines)
        for part in lines:
            print("part: ", part)
            if part.strip():
                parts = part.split(":")
                key = parts[0].strip()
                value = parts[1].strip()
                results.append({key: value})

        return results


def run_2019_learna(target_structure=None, target_id=1, timeout=600, learning_rate=0.0005991629320464973,
                    rl_algo_type="bin/learna", input_data_dir=None, hamming_tol=10,
                    no_shared_agent=False):

    """
    run_learna
    ----------
    Puts the parameters needed to run the LEARNA algorithms into a subprocess.
    This allows the function to be easily called through Python, and the output
    can be utilised in a Python script.

    rl_algo_type: Can be one of the four following strings:
                    bin/learna: Runs the LEARNA algorithm
                    bin/meta-learna: Runs the Meta-LEARNA algorithm
                    bin/meta-learna-adapt: Runs the Meta-LEARNA-Adapt algorithm
    """
    error = False
    learna_script = rl_algo_type
    if (target_structure == None and input_data_dir == None):
        print("You must input either the target structure, or the"
              "pathway to the target structures")
    if rl_algo_type=="bin/liblearna":
        args = [
            "--input_file", input_data_dir,
            "--timeout", str(timeout),
            "--learning_rate", str(learning_rate),
            "--hamming_tolerance", str(hamming_tol),
        ]
    elif target_structure != None:

        args = [
            "--target_structure", str(target_structure),
            "--timeout", str(timeout),
            "--learning_rate", str(learning_rate),
            # "--data_dir", str(data_dir),
            # "--dataset", str("eterna")
            "--hamming_tolerance", str(hamming_tol),
            "--mutation_threshold", str(5),
            "--reward_exponent", str(8.93),
            "--state_radius", str(29),
            "--conv_sizes", str([11, 3]),
            "--conv_channels", str([10, 3]),
            "--num_fc_layers", str(1),
            "--fc_units", str(52),
            "--batch_size", str(123),
            "--entropy_regularization", str(1.51e-4),
            "--lstm_units", str(3),
            "--num_lstm_layers", str(0),
            "--embedding_size", str(2),
            "--batch_size", str(126),
        ]
    else:
        args = [
            "--timeout", str(timeout),
            "--learning_rate", str(learning_rate),
            "--input_file", input_data_dir,
            "--hamming_tolerance", str(hamming_tol),
            # "--dataset", str("eterna")
            # TODO: Add more arguments
        ]

    # Construct the command to execute learna.py with arguments
    command = ["python", learna_script] + args

    print("Running Learna Command...")
    print(command)
    try:
        output = subprocess.check_output(command)
        # Process the output
        output_string = output.decode().replace("\\n", "\n")
        # print("Output:", output_string)

    except subprocess.CalledProcessError as e:
        output_string = e.output.decode().replace("\\n", "\n")
        print("Error occurred while executing:", output_string)
        error = True

    if error:
        return output_string
    results = []
    lines = output_string.split("\n")
    for line in lines:
        if "### Number of targets:" in line:
            targets = int(line.split(":")[1].strip())
            break

    for i in range(1, targets+1):
        # Find the index of the line containing the headers
        try:
            header_index = lines.index("Solutions for target structure "+str(i))
        except:
            new_row = {"Sequence": "No match found", "Hamming Distance": 1}
            results.append(new_row)
            continue


        data_lines = lines[header_index + 3:]  # Skip the header and separator lines
        # Extract sequence and hamming_distance values for each row
        x = 0
        for line in data_lines:
            if line.__contains__("-"):
                continue
            elif line.strip():  # Skip empty lines
                # print("LINE: ", line)
                parts = line.split('|')
                sequence = parts[6].strip()
                hamming_distance = float(parts[4].strip())
                rel_hamming_distance = float(parts[5].strip())
                new_row = {"Sequence": sequence, "Hamming Distance": hamming_distance, \
                           "Normalised Hamming Distance": rel_hamming_distance}
                results.append(new_row)
                break

    return results

def make_model(input_data_dir, timeout, worker_count, save_path,
               restore_path, learning_rate=5.99e-4,
               reward_exponent=9.34, state_radius=32, conv_sizes=[17, 5],
               conv_channels=[7, 18], fc_layers=1, fc_units=57,
               batch_size=126, entropy_reg=6.76e-5, embed_size=3,
               lstm_units=1, lstm_layers=28):
    args = [
        "--timeout", str(timeout),
        "--data_dir", input_data_dir,
        "--worker_count", str(worker_count),
        "--save_path", save_path,
        "--restore_path", restore_path,
        "--dataset", str("rfam_learn_train"),
        # "--target_structure_ids", str(target_structure_ids),
        # TODO: Add more arguments
    ]

    command = ["python", "learna_tools/learna/learn_to_design_rna.py"] + args

    print("Running Learna Command...")
    print(command)
    try:
        output = subprocess.check_output(command)
        # Process the output
        output_string = output.decode().replace("\\n", "\n")
        # print("Output:", output_string)

    except subprocess.CalledProcessError as e:
        output_string = e.output.decode().replace("\\n", "\n")
        print("Error occurred while executing:", output_string)
        error = True


def add_to_seq_ratio(sequence, the_dict):
    print(the_dict)
    if len(sequence) >= 350:
        the_dict["350+"] += 1
    elif len(sequence) >= 150:
        the_dict["150-349"] += 1
    elif len(sequence) >= 100:
        the_dict["100-149"] += 1
    elif len(sequence) >= 50:
        the_dict["50-99"] += 1
    else:
        the_dict["0-49"] += 1
    return the_dict

if __name__ == '__main__':

    target_structure = "...(((((....)))))..."
    data_directory = "Eterna100.txt"
    # liblearna_data_dir = "/Users/jackcole/PycharmProjects/learna_tools2/learna_tools/liblearna/data/data/eterna/raw/eterna_raw.txt"
    # target_structure_ids = 1
    # make_model(input_data_dir=data_directory, timeout=10, worker_count=4,
    #            save_path="models/", restore_path="models/")

    # result = run_2019_learna(target_structure=target_structure, rl_algo_type="bin/learna", timeout=10)
    """learna_result = run_2019_learna(input_data_dir=data_directory,
                                    rl_algo_type="bin/learna", timeout=3600,
                                    hamming_tol=0)
    print(learna_result)
    """
    meta_learna_result = run_2019_learna(input_data_dir=data_directory,
                                         rl_algo_type="bin/meta-learna",
                                         timeout=3600,hamming_tol=0,
                                         learning_rate=6.44e-5)
    """
    meta_learna_adapt_result = run_2019_learna(input_data_dir=data_directory,
                                               rl_algo_type="bin/meta-learna-adapt",
                                               timeout=3600, hamming_tol=1,
                                               learning_rate=6.44e-5)
    
    lib_learna_result = run_2019_learna(input_data_dir=liblearna_data_dir,
                                               rl_algo_type="bin/liblearna",
                                               timeout=100, hamming_tol=1,
                                               learning_rate=5.9e-4)
    """
    #mcts_result = run_mcts(input_data_dir=data_directory, timeout=600)
    #df = pd.read_csv(data_directory)
    # sequences = df[["str"]]
    # print(sequences)

    # print("LEARNA RESULT: ", learna_result)
    learna_number_of_success = 0
    learna_avg_seq_len = 0
    learna_seq_ratios = {"0-49": 0, "50-99": 0, "100-149": 0, "150-349": 0, "350+": 0}
    meta_learna_seq_ratios = {"0-49": 0, "50-99": 0, "100-149": 0, "150-349": 0, "350+": 0}
    meta_learna_adapt_seq_ratios = {"0-49": 0, "50-99": 0, "100-149": 0, "150-349": 0, "350+": 0}
    mcts_seq_ratios = {"0-49": 0, "50-99": 0, "100-149": 0, "150-349": 0, "350+": 0}
    meta_learna_number_of_success = 0
    meta_learna_avg_seq_len = 0
    meta_learna_adapt_number_of_success = 0
    meta_learna_adapt_avg_seq_len = 0
    mcts_number_of_success = 0
    mcts_seq_len = 0
    """
    for rna_dict in learna_result:
        if rna_dict["Hamming Distance"] == 0:
            learna_number_of_success += 1
            learna_avg_seq_len += len(rna_dict["Sequence"])
            learna_seq_ratios = get_seq_ratio(rna_dict["Sequence"], learna_seq_ratios)
    #learna_avg_seq_len = learna_avg_seq_len / learna_number_of_success
    """
    for rna_dict in meta_learna_result:
        print(rna_dict)
        if rna_dict["Hamming Distance"] == 0:
            meta_learna_number_of_success += 1
            meta_learna_avg_seq_len += len(rna_dict["Sequence"])
            meta_learna_seq_ratios = add_to_seq_ratio(rna_dict["Sequence"], meta_learna_seq_ratios)
    meta_learna_avg_seq_len = meta_learna_avg_seq_len / meta_learna_number_of_success
    """
    for rna_dict in meta_learna_adapt_result:
        print(rna_dict)
        if rna_dict["Hamming Distance"] == 0:
            meta_learna_adapt_number_of_success += 1
            meta_learna_adapt_avg_seq_len += len(rna_dict["Sequence"])
            meta_learna_adapt_seq_ratios = get_seq_ratio(rna_dict["Sequence"], meta_learna_adapt_seq_ratios)
    meta_learna_adapt_avg_seq_len = meta_learna_adapt_avg_seq_len / meta_learna_adapt_number_of_success
    for x in range(len(mcts_result)):
        print(mcts_result[x])
        print(mcts_result[x]["Normalised Hamming Distance"])
        if mcts_result[x]["Normalised Hamming Distance"] == "1.0":
            mcts_number_of_success += 1
            mcts_seq_len += len(mcts_result[x]["Sequence"])
            mcts_seq_ratios = get_seq_ratio(mcts_result[x]["Sequence"], mcts_seq_ratios)
    mcts_seq_len = mcts_seq_len / mcts_number_of_success
    
    print("LEARNA SUCCESSES: ", learna_number_of_success)
    print("LEARNA AVG SEQ LENGTH: ", learna_avg_seq_len)
    print("LEARNA SEQ UNDER 50 SOLVED: ", learna_seq_ratios["0-49"])
    print("LEARNA SEQ UNDER 100 SOLVED: ", learna_seq_ratios["50-99"])
    print("LEARNA SEQ UNDER 150 SOLVED: ", learna_seq_ratios["100-149"])
    print("LEARNA SEQ ABOVE 150 SOLVED: ", learna_seq_ratios["150+"])
    """
    print("META-LEARNA SUCCESSES: ", meta_learna_number_of_success)
    print("META-LEARNA AVG SEQ LENGTH: ", meta_learna_avg_seq_len)
    print("META-LEARNA SEQ UNDER 50 SOLVED: ", meta_learna_seq_ratios["0-49"])
    print("META-LEARNA SEQ UNDER 100 SOLVED: ", meta_learna_seq_ratios["50-99"])
    print("META-LEARNA SEQ UNDER 150 SOLVED: ", meta_learna_seq_ratios["100-149"])
    # print("META-LEARNA SEQ ABOVE 150 SOLVED: ", meta_learna_seq_ratios["150+"])
    df = pd.DataFrame(meta_learna_result)
    df.to_csv("M_Learna_Result_1_Hour.csv")

    """
    print("META-LEARNA-ADAPT SUCCESSES: ", meta_learna_adapt_number_of_success)
    print("META-LEARNA-ADAPT AVG SEQ LENGTH: ", meta_learna_adapt_avg_seq_len)
    print("META-LEARNA-ADAPT SEQ UNDER 50 SOLVED: ", meta_learna_adapt_seq_ratios["0-49"])
    print("META-LEARNA-ADAPT SEQ UNDER 100 SOLVED: ", meta_learna_adapt_seq_ratios["50-99"])
    print("META-LEARNA-ADAPT SEQ UNDER 150 SOLVED: ", meta_learna_adapt_seq_ratios["100-149"])
    print("META-LEARNA-ADAPT SEQ ABOVE 150 SOLVED: ", meta_learna_adapt_seq_ratios["150+"])
    print("MCTS SUCCESSES: ", mcts_number_of_success)
    print("MCTS AVG SEQ LENGTH: ", mcts_seq_len)
    print("MCTS SEQ UNDER 50 SOLVED: ", mcts_seq_ratios["0-49"])
    print("MCTS SEQ UNDER 100 SOLVED: ", mcts_seq_ratios["50-99"])
    print("MCTS SEQ UNDER 150 SOLVED: ", mcts_seq_ratios["100-149"])
    print("MCTS SEQ ABOVE 150 SOLVED: ", mcts_seq_ratios["150+"])
"""
"""

    # To use LibLEARNA we require more information within the datasets:
    # id, structure constraints, sequence constraints, GC content, desired energy


    # NOTE: THERE WAS AN UNDERSCORE IN THE 77TH ETERNA100 SEQUENCE MAKING IT CRASH
    # I HAVE REMOVED THIS UNDERSCORE TO PREVENT IT FROM BREAKING - NOT SURE IF IT
    # REPRESENTS A PSEUDOKNOT OR IS JUST A MISTAKE

"""