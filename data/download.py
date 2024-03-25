import gdownloader as gdown
import subprocess
import pandas as pd

from pathlib import Path


DATASETS = {
  'ArchiveII_pk.plk.gz': 'https://drive.google.com/uc?id=1W8hbU03lCdttuh5zKXPlU5aCHdk8p2aY',  # https://drive.google.com/file/d/1W8hbU03lCdttuh5zKXPlU5aCHdk8p2aY/view?usp=drive_link
  'anta_pseudo-test.plk.gz': 'https://drive.google.com/uc?id=1ORh9HUJkZvMg2UoufvrEoXOpn0JceUS_',  # https://drive.google.com/file/d/1ORh9HUJkZvMg2UoufvrEoXOpn0JceUS_/view?usp=drive_link
  'anta_rfam-test.plk.gz': 'https://drive.google.com/uc?id=1Fpl1reKaLejgVbJbMxqrmDDj7NsWF3zD',  # https://drive.google.com/file/d/1Fpl1reKaLejgVbJbMxqrmDDj7NsWF3zD/view?usp=drive_link
  'anta_strand-test.plk.gz': 'https://drive.google.com/uc?id=1JnRJIci07h_OptZGgwkzu8HsCxQw567G',  # https://drive.google.com/file/d/1JnRJIci07h_OptZGgwkzu8HsCxQw567G/view?usp=drive_link
  'rfam_learn-test.plk.gz': 'https://drive.google.com/uc?id=1JlFIfljDd8pvNfKkPyX8di9Qzpbp-GvB',  # https://drive.google.com/file/d/1JlFIfljDd8pvNfKkPyX8di9Qzpbp-GvB/view?usp=drive_link
  'rfam_taneda-test.plk.gz': 'https://drive.google.com/uc?id=1D1_0uGM-z1pASHKtdnTVqj8fzCxlsh72',  # https://drive.google.com/file/d/1D1_0uGM-z1pASHKtdnTVqj8fzCxlsh72/view?usp=drive_link
  'eterna100_v1.plk.gz': 'https://drive.google.com/uc?id=1ciA2vIFZ10ukAeZ88-sRqpTHyfutfgjD',  # https://drive.google.com/file/d/1ciA2vIFZ10ukAeZ88-sRqpTHyfutfgjD/view?usp=drive_link
  'eterna100_v2.plk.gz': 'https://drive.google.com/uc?id=1Iry9e_VFsE_hXLAlIuycV1VpmfnF49rJ',  # https://drive.google.com/file/d/1Iry9e_VFsE_hXLAlIuycV1VpmfnF49rJ/view?usp=drive_link
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
    # Change to "wget" rather than "usr/local/" if it doesn't work
    subprocess.call(["/usr/local/bin/wget", download_link], cwd=destination)
    # subprocess.call(["unzip", file_id], cwd=destination)
    # subprocess.call(["rm", file_id], cwd=destination)

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
        #print('### Download Rfam family CMs')
    #if 'CMs' in args.data_types:
    #    Path(args.rfam_cm_outdir).mkdir(exist_ok=True, parents=True)
    #    download_rfam_cms(args.rfam_cm_outdir)
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
