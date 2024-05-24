import subprocess

def fold_with_pKiss(seq):
    # Construct the command
    cmd = ["pKiss", "--strategy", "P", "--mode", "mfe", seq]

    try:
        # Run the command and capture the output
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        output = result.stdout.decode().strip().split()
        print(output)
        if (len(output) > 1):
            mfe = output[-2]
            structure = output[-1]
            return [structure, mfe], [0] * len(output)
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Error occurred: {e}")
        return None


def ppKiss(listOfSeqs):
    result = []
    for s in listOfSeqs:
        result.append(fold_with_pKiss(s))
    result = array(result)

    return result, [0] * len(result)

if __name__ == '__main__':
    seq = "AUAUGCGCGCGCAUAU"
    print(fold_with_pKiss(seq))
