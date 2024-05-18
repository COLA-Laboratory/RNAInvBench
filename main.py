from tqdm import tqdm
import threading

# algorithms
from rna_design_algorithms.random_based import rnassd
from rna_design_algorithms.sample_based.IncaRNAtion import incarnation
from rna_design_algorithms import rnainverse

# evaluations
from evaluations.successRate import successRate
from evaluations.structureDistance import structureDistance
from evaluations.MFEgcContent import MFE_GCcontent

# datasets
from utilis.select_dataset import select_dataset
from utilis.check_and_mk import check_and_mk

def IncaRNAtion(targets, folder, dataset):
    pbar = tqdm([i for i in range(len(targets))])
    for i in pbar:
        target = targets[i]
        pbar.set_description(f'IncaRNAtion: Processing structure {i}, length: {len(target)}')
        incarnation_save = folder + "incarnation_" + dataset + "_" + str(i)
        incarnation.call_incarnation(target, save_file=incarnation_save)

def RNA_SSD(targets, folder, dataset):
    pbar = tqdm([i for i in range(len(targets))])
    for i in pbar:
        target = targets[i]
        pbar.set_description(f'RNA-SSD: Processing structure {i}, length: {len(target)}')
        rnassd_save = folder + "rnassd_" + dataset + "_" + str(i)
        rnassd.call_rnassd(target, save_file=rnassd_save)

def RNAinverse(targets, folder, dataset):
    pbar = tqdm([i for i in range(len(targets))])
    for i in pbar:
        target = targets[i]
        pbar.set_description(f'RNAinverse: Processing structure {i}, length: {len(target)}')
        rnainverse_save = folder +"rnainverse_" + dataset + "_" + str(i)
        rnainverse.call_rnainverse(target, save_file=rnainverse_save)

if __name__ == "__main__":
    # choose database
    dir = 'data/'
    file = 'inverse_rna_folding_benchmark_dotbracket.pkl.gz'
    path = dir + file
    dataset = 'eterna100_v2'
    range_len = [1, 500]

    targets, _ = select_dataset(path, dataset, range_len)

    folder = 'results/' + dataset + '/'
    check_and_mk(folder)

    # design
    # 1. incarnation
    thread_incarnation = threading.Thread(target=IncaRNAtion, args=(targets, folder, dataset))
    thread_rnainverse = threading.Thread(target=RNAinverse, args=(targets, folder, dataset))

    # start threads
    thread_incarnation.start()
    thread_rnainverse.start()

    # completed.
    thread_incarnation.join()
    thread_rnainverse.join()
    print("all mission completed.")
    
    # # evaluate
    # res = successRate(60, result_folder=folder)
    # for key, value in res.items():
    #     print(f'{key}: {value}')

    # structure_df = structureDistance(plot=True, result_folder=folder)

    # df_MFEgc = MFE_GCcontent(result_folder=folder)


