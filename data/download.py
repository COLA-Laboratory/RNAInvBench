import gdownloader as gdown
import subprocess
import wget
import pandas as pd

from pathlib import Path


DATASETS = {
  'constrained_design_benchmark.plk.gz': 'https://drive.google.com/uc?id=1ZU_VKBveT25GMuEyMxK6_F9DcMp_N1bY',  # https://drive.google.com/file/d/1I9o3QDHWYG3GjrqzPlpw7h_nWA-ld2Zg/view?usp=drive_link; https://drive.google.com/file/d/1ZU_VKBveT25GMuEyMxK6_F9DcMp_N1bY/view?usp=drive_link
  'constrained_design_train.plk.gz': 'https://drive.google.com/uc?id=12irVqG4al3UaYdfEQArQoKOChUhJGGus',  # https://drive.google.com/file/d/1UDZ8uHpbcRXZla2opQFa0D8eOuQs_9-i/view?usp=drive_link; https://drive.google.com/file/d/1NsmmviFcJ86EdsC9eE5TVL6Nm5nS5oks/view?usp=drive_link; https://drive.google.com/file/d/12irVqG4al3UaYdfEQArQoKOChUhJGGus/view?usp=drive_link
  'constrained_design_valid.plk.gz': 'https://drive.google.com/uc?id=1MeXC2E5Norf2rT4yXbeSo9xQB4NP2IrN',  # https://drive.google.com/file/d/1FXksh4E789mjQWYOmhoRN-nxt_km1j72/view?usp=drive_link; https://drive.google.com/file/d/1MeXC2E5Norf2rT4yXbeSo9xQB4NP2IrN/view?usp=drive_link
  'inverse_rna_folding_benchmark.plk.gz': 'https://drive.google.com/uc?id=1mSH_upJHBlC9XoUcMVH1MucyVLt4Rhgx',  # https://drive.google.com/file/d/1hlcKaPLjgfdZB2SfvjEmCO6qcK3DC65e/view?usp=drive_link; https://drive.google.com/file/d/1mSH_upJHBlC9XoUcMVH1MucyVLt4Rhgx/view?usp=drive_link
  'inverse_rna_folding_train.plk.gz': 'https://drive.google.com/uc?id=1baVwuEWJnBlFJWrypWE6AgPJzP2reUYw',  # https://drive.google.com/file/d/1dj1SWCJ3bxHge__yrD6VH2apGAHAOWff/view?usp=drive_link; https://drive.google.com/file/d/1_BIrEx19MQINIB2Ge3LabUyJMqcR_k6H/view?usp=drive_link; https://drive.google.com/file/d/1baVwuEWJnBlFJWrypWE6AgPJzP2reUYw/view?usp=drive_link
  'inverse_rna_folding_valid.plk.gz': 'https://drive.google.com/uc?id=1-3UBqU4WPvH_hEJTaekx_U7QAMcSLKW9',  # https://drive.google.com/file/d/1A6Tsoh8mJM9c1gho8_h5wOesbWVpguVZ/view?usp=drive_link; https://drive.google.com/file/d/1-3UBqU4WPvH_hEJTaekx_U7QAMcSLKW9/view?usp=drive_link
}

def select_and_download(task, save_dir, all_data=True):
    if all_data:
        for s in ['_train', '_valid', '_benchmark']:
            try:
                download_single_file_from_gdrive(url=DATASETS[f"{task}{s}.plk.gz"],
                                                 destination=save_dir + '/' + f"{task}{s}.plk.gz")
            except Exception as e:
                pass
    else:
        try:
            download_single_file_from_gdrive(url=DATASETS[f"{task}_benchmark.plk.gz"],
                                             destination=save_dir + '/' + f"{task}_benchmark.plk.gz")
        except Exception as e:
            pass


def download_single_file_from_gdrive(url, destination):
    gdown.download(url, destination, quiet=False)

def download_3d_data(
                     destination : str,
                     resolution : str = '3.5',  # available options: 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 20.0, all
                     data_type : str = 'cif',  # available options: cif, fasta, pdb
                     bgsu_version : float = '3.286',
                     rna_type : str = 'solo',  # available options: 'solo' -> RNA only data; 'all' -> all molecules; 'complex' -> RNA-protein complexes; 'hybrid' -> DNA-RNA hybrids
                     representatives_only : bool = False,  # if True: download only representative members of every RNA equivalence class
                     ):
    rep = 'representative' if representatives_only  else "member"
    bgsu = str(bgsu_version).replace('.', '_')
    resolution = resolution.replace('.', '_')

    download_link = f'https://rnasolo.cs.put.poznan.pl/media/files/zipped/bunches/{data_type}/{rna_type}_{rep}_{data_type}_{resolution}__{bgsu}.zip'
    file_id = f"{rna_type}_{rep}_{data_type}_{resolution}__{bgsu}.zip"
    destination = f"{destination}/{rna_type}_{rep}_{data_type}_{resolution}__{bgsu}"
    Path(destination).mkdir(exist_ok=True, parents=True)
    # Change to ["wget"] if this doesn't work - current solution works on MAC
    subprocess.call(["/usr/local/bin/wget", download_link], cwd=destination)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_types', type=str, default=['2D', 'rnasolo'], nargs="+", help='Define which data to download')  # also available: 'rnasolo'
    parser.add_argument('--bgsu', default=3.326, type=float, help='download all 3D rna data from RNAsolo')
    parser.add_argument('--data_outdir', default='data', type=str, help='output_directory for 2D data')
    parser.add_argument('--rfam_cm_outdir', default='RnaBench/lib/data/CMs', type=str, help='output_directory for Rfam.cm')

    args = parser.parse_args()
    if '2D' in args.data_types:
        print('### Downloading all secondary structure data.')
        Path(args.data_outdir).mkdir(exist_ok=True, parents=True)
        for k, v in DATASETS.items():
            destination = f'{args.data_outdir}/{k}'
            if not Path(destination).is_file():
                url = v
                download_single_file_from_gdrive(url, destination)
            else:
                print(f'File {destination} already exists. Skipping.')
    if 'rnasolo' in args.data_types:
        print('### Download 3D data')
        bgsu = args.bgsu
        threedee_destination = Path(args.data_outdir, '3D')
        threedee_destination.mkdir(exist_ok=True, parents=True)
        res = ['1.5', '2.0', '2.5', '3.0', '3.5', '4.0', '20.0', 'all']
        rep = [True, False]
        d_list = [[r, rep[0]] for r in res] + [[r, rep[1]] for r in res]
        for r, rep in d_list:
            try:
                download_3d_data(str(threedee_destination.resolve()),
                representatives_only=rep,
                resolution=r,
                bgsu_version=bgsu,
                )
            except Exception as e:
                print(e)
                print('### 3D data download not working. Try to increase the BGSU version via --bgsu.', f"Current version is {bgsu}")
                continue
