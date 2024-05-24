from grakel import Graph
from grakel.kernels import (
    WeisfeilerLehman,
    VertexHistogram,
    #WeisfeilerLehmanOptimalAssignment,
    #ShortestPath,
    #Propagation,
    #RandomWalkLabeled,
    #CoreFramework,
)

from learna_tools.metrics.matrix_metrics import db2mat


def weisfeiler_lehmann(designed_sequence=None, designed_structure=None, target_sequence=None, target_structure=None):

    n_graph_iter = 5
    wl_kernel = WeisfeilerLehman(n_iter=n_graph_iter, normalize=True, base_graph_kernel=VertexHistogram)
    true_mat = db2mat(target_structure)
    pred_mat = db2mat(designed_structure)
    true_graph = Graph(initialization_object=true_mat.astype(int), node_labels={s:str(s) for s in range(true_mat.shape[0])})
    pred_graph = Graph(initialization_object=pred_mat.astype(int), node_labels={s:str(s) for s in range(pred_mat.shape[0])})

    wl_kernel.fit_transform([true_graph])
    wl_pred = wl_kernel.transform([pred_graph])

    return 1 - wl_pred[0][0]
