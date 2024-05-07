import re

import requests
from bs4 import BeautifulSoup

def call_rnassd(target_structure, constraints='', temperature=37, tries=1, GCcontent_paired=50, GCcontent_unpaired=50, seed=''):
    '''
    Call RNA-SSD from RNA Designer online web application to find an RNA sequence that folds into the target structure.
    Using ViennaRNA packages inside on secondary structure prediction.

    Parameters:
    - target_structure (str): The target RNA secondary structure in dot-bracket notation.
        - maximum length: 500
    - seq_constraints (str): Sequence Constraints of targeted strcture.
        - maximum length: 500
    - temperature (float): temperature of ViennaRNA secondary structure prediction.
        - range: [0.0, 100.0]
    - tries (int): The number of sequences that are to be designed.
        - max: 10
    - GCcontent_paired/GCcontent_unpaired(float): be used to control the target GC content of paired/unpaired bases in the desired structure, and are specified as percentages
        - range: [0, 100]
    - seed (int): random seed

    Returns:
    - output_nSeq: a list consists of designed RNA outputs
        dict keys:
        - designed_sequence (str)
        - designed_structure (str)
        - distance (int): Distance to desired structure
        - MFE (float): Free energy of the predicted structure in kcal/mol, as computed by the Vienna RNA Package RNAfold algorithm
        - MFE_ther (float): Free energy of the thermodynamic ensemble in kcal/mol
        - MFE_prob (float): Probability of the minimal free energy structure in the thermodynamic ensemble
        - gc_paired/gc_unpaired (float): GC content of paired/unpaired regions
        - seed (int): random seed used in this design loop
    '''
    # post
    web_url = 'http://www.rnasoft.ca/cgi-bin/RNAsoft/RNAdesigner/rnadesign.pl'

    data = {
        'structure': target_structure,
        'constraints': constraints,
        'temperature': str(temperature),
        'tries': str(tries),
        'GCcontent_paired': str(GCcontent_paired),
        'GCcontent_unpaired': str(GCcontent_unpaired),
        'seed': str(seed),
        'bywebpage': 'checked',
        'Submit': 'Run RNA Designer'
    }

    # response
    response = requests.post(web_url, data=data)

    # extract results
    soup = BeautifulSoup(response.text, 'html.parser')
    terms = soup.select('tr')

    # input parameters
    input_dict = {}
    input_trs = terms[3:10]
    for tr in input_trs:
        key = tr.font.text.split('[')[0].strip()
        value = tr.b.text.strip()

        if key == 'Structure':
            input_dict['target_structure'] = value
        elif key == 'Sequence constraints':
            input_dict['seq_constrints'] = value
        elif key == 'Temperature':
            input_dict['temperature'] = float(value.split('°')[0])
        elif key == 'Number of sequences to design':
            input_dict['n_seq'] = int(value)
        elif key == 'Target GC content in paired regions':
            input_dict['gc_paired'] = float(value.split('%')[0])
        elif key == 'Target GC content in unpaired regions':
            input_dict['gc_unpaired'] = float(value.split('%')[0])
        elif key == 'Random number seed':
            input_dict['seed'] = int(value)

    # output parameters
    output_nSeq = []
    start_ind = 11
    len = 8
    for i in range(input_dict['n_seq']):
        # for every tries
        output_dict = {}
        output_trs = terms[start_ind: start_ind + len]

        for tr in output_trs:
            key = tr.font.text.split('[')[0].strip()
            value = tr.b.text.strip()

            if re.search(r'Sequence \d+\s*and MFE structure', key):
                designed_sequence = re.search(r"([A-Z]+)", value).group()
                designed_structure = re.search(r"([()\.]+)", value).group()
                output_dict['designed_sequence'] = designed_sequence
                output_dict['designed_structure'] = designed_structure
            elif key == 'Distance to desired structure':
                output_dict['distance'] = int(value)
            elif key == 'Minimum free energy':
                output_dict['MFE'] = float(value.split(' ')[0])
            elif key == 'Free energy of the thermodynamic ensemble':
                output_dict['MFE_ther'] = float(value.split(' ')[0])
            elif key == 'Probability of MFE in ensemble':
                output_dict['MFE_prob'] = float(value)
            elif key == 'Actual GC content of paired regions':
                output_dict['gc_paired'] = float(value.split('%')[0])
            elif key == 'Actual GC content of unpaired regions':
                output_dict['gc_unpaired'] = float(value.split('%')[0])
            elif key == 'Random number seed for this sequence':
                output_dict['seed'] = int(value)
        # add to outputs list
        output_nSeq.append(output_dict)
        start_ind += 1

    return output_nSeq

# example usage
if __name__ == "__main__":
    target_structure = "(((....)))((...))"
    outputs = call_rnassd(target_structure, tries=3)
    n = 0  # may have many designed sequences(here is 3)
    for i in outputs:
        n += 1
        print(f'Sequence {n}(order by MFE)')
        for key, value in i.items():
            print(f'{key}: {value}')
        print('----------------')

