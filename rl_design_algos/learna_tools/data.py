import pandas as pd

from pathlib import Path

class InputHandler():
    def __init__(self, input_dir=None, input_file=None, format=None):
        pass


class OutputHandler():
    def __init__(self, out_dir=None, out_file=None, format=None):
        pass

def read_fasta(path):
    pass

def write_fasta(data, outpath):
    pass

def read_task_description(path):
    """
    Reads task description from Path.

    tasks descriptions files have the following format:

    RNAid
    sequence constraints
    structure constraints
    ?local gc constraints?  # add later...
    global gc content
    global energy

    Example:
    > Test RNA
    #seq ACGCGGCGCUA AGUUU UUCGCA NNNNNNN  # sequence conststraints
    #str ..(((...))) NNN.. ||)..( NNNNNN|  # structure constraints
    #localgc - - - 0.2                         # partial gc-contents
    #globgc 0.6                               # global gc-content
    #globen -13.4                             # global energy
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    tasks = []
    set_as_df = False
    for i, line in enumerate(lines):
        if line.startswith('>'):
            line = line.replace('\n', '')
            task = {'Id': line[1:]}
            for l in lines[i+1:]:
                if not l.strip():
                    continue

                # l = l.rstrip()
                if l.startswith('>'):
                    break
                if l.startswith('#'):
                    l_list = l.split()
                    if l_list[0][1:] == 'str':
                        l = ''.join(l.split('#str ')).replace('\n', '')
                        task['str'] = l.replace(' ', 'X')
                    elif l_list[0][1:] == 'seq':
                        l = ''.join(l.split('#seq ')).replace('\n', '')
                        task['seq'] = ''.join(l_list[1:]).replace('X', '').replace(' ', '')
                    else:
                        task[l_list[0][1:]] = ''.join(l_list[1:])
                else:
                    task['str'] = l.rstrip()
            tasks.append(task)
        else:
            set_as_df = True

    if set_as_df:
        data = pd.read_csv(path)
    else:
        data = pd.DataFrame(tasks)

    return data



def read_fastas_from_dir():
    pass


if __name__ == '__main__':
    tmp_path = 'test_file.in'

    with open(tmp_path, 'w+') as f:
        for i in range(20):
            f.write(f"> {i}\n")
            f.write("#seq ACGCGGCGCUA AGUUU UUCGCA NNNNNNN\n")
            f.write("#str ..(((...))) NNN.. ||)..( NNNNNN|\n")
            f.write("#globgc 0.6\n")
            f.write("#globen -13.4\n")
            if i % 2 == 0:
                f.write("#localgc - - - 0.2\n")

    tasks = read_task_description(tmp_path)
    print(tasks)
    Path(tmp_path).unlink()

    with open(tmp_path, 'w+') as f:
        for i in range(20):
            f.write(f"> {i}\n")
            f.write("...((((...))))...\n")


    tasks = read_task_description(tmp_path)
    print(tasks)
    Path(tmp_path).unlink()