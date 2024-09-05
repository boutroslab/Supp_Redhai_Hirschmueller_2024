import pandas as pd
import numpy as np
import graphtools as gt
import phate
import scprep
import meld
import scipy


pca_data = pd.read_csv("./output/MELD/pca_data_tx23.tsv", sep="\t",
                       index_col="cell_barcode")
mdata = pd.read_csv(".output/MELD/mdata_tx23.tsv", sep="\t",
                    index_col="cell_barcode", low_memory=False)
                    
data = pd.read_csv("./output/MELD/norm_cnts_tx23.tsv", sep="\t",
                   index_col="cell_barcode")
data = scipy.sparse.csr_matrix(data)

benchmarker = meld.Benchmarker()

benchmarker.fit_phate(pca_data)

from joblib import Parallel, delayed


def simulate_pdf_calculate_likelihood(benchmarker, seed, beta):
    benchmarker.set_seed(seed)
    benchmarker.generate_ground_truth_pdf()

    benchmarker.generate_sample_labels()
    benchmarker.calculate_MELD_likelihood(beta=beta)
    MELD_mse = benchmarker.calculate_mse(benchmarker.expt_likelihood)
    return MELD_mse, seed, beta, benchmarker.graph.knn


knn_range = np.arange(3, 20)
beta_range = np.arange(10, 150)

results = []

with Parallel(n_jobs=65) as p:
    for knn in knn_range:
        # doing this outside the parallel loop because building the graph takes the longest
        benchmarker.fit_graph(data, knn=knn)
        print(knn)
        curr_results = p(delayed(simulate_pdf_calculate_likelihood)(benchmarker, seed, beta) \
                         for seed in range(2) for beta in beta_range)
        curr_results = pd.DataFrame(curr_results, columns=['MSE', 'seed', 'beta', 'knn'])
        results.append(curr_results)

results = pd.concat(results, axis=0)
results.to_csv("./output/MELD/parameter_search_meld_tx23.csv", sep="\t",index=False)
