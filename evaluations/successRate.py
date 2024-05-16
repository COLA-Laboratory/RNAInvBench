import pandas as pd
import os

def successRate(time_limit, result_folder='results'):
    """
    success rate in limited time.

    parameters:
    - result_folder (str): folder path of results of one algorithm.
    - time_limit (float): time limitation of execution

    return:
    rate (float): success rate
    """
    successRate_dict = {}
    file_list = os.listdir(result_folder)

    # group file by method name
    method_groups = {}
    for file_name in file_list:
        method_name = file_name.split('_')[0]
        if method_name not in method_groups:
            method_groups[method_name] = []
        method_groups[method_name].append(file_name)

    for method, files in method_groups.items():
        n_total = len(files)
        n_success = 0
        for file in files:
            if file.endswith('.pkl'):
                file_path = os.path.join(result_folder, file)
                df = pd.read_pickle(file_path)  # read the result pickle
                try:
                    if df['sequence'][0] and float(df['time'][0]) <= time_limit:
                        n_success = n_success + 1
                except KeyError:
                    continue
        # rate
        successRate_dict[method] = n_success / n_total

    return successRate_dict



