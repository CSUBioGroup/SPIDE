import pandas as pd
import numpy as np
from tqdm import tqdm

# from rpy2.robjects import r
# import rpy2.robjects as robjects

from numpy import linalg as LA
from utils import *
from algo import *

import concurrent.futures
from itertools import repeat


def spide(geneexp_file, k, ncore, ppi_file, result_file):
    ppi_matrix, gene_id, id_gene = extract_ppi(ppi_file)

    print("PPI Loaded")

    expdf = pd.read_csv(geneexp_file, header=None)

    gene_exp = {}
    for idx, row in expdf.iterrows():
        gname = int(row[0])
        exps = row[1:].tolist()
        gene_exp[gname] = exps

    print("GeneExp Loaded")

    gene_exp, ppi_matrix, gene_id, id_gene = DoIntegPPI(gene_exp, ppi_matrix, gene_id, id_gene)
    lgene_exp = np.log2(gene_exp+1.1)

    print("Solving Eig")
    w,v = LA.eig(ppi_matrix-np.identity(ppi_matrix.shape[0]))
    maxSR = np.log(max(w).real)
    print("Eig Solved")

    cellpear = np.corrcoef(lgene_exp.T)
    row, col = np.diag_indices_from(cellpear)
    cellpear[row,col] = 0


    fexp = fillfunc2(gene_exp, cellpear, k)

    tempdf = pd.DataFrame(fexp)
    tempdf = quantileNormalize(tempdf)
    fexp = tempdf.values
    fexp = np.log2(fexp+1.1)

    print("ProcessPool")
    with concurrent.futures.ProcessPoolExecutor(ncore) as p:
        cells = range(cellpear.shape[0])
        ans = p.map(newSCENT2ps_abs_fill_sub, repeat(fexp), repeat(ppi_matrix), repeat(cellpear), repeat(k), repeat(lgene_exp), cells, repeat(maxSR))
    # p = concurrent.futures.ProcessPoolExecutor(10)
    # cells = range(cellpear.shape[0])
    # ans = p.map(newSCENT2ps_abs_fill_sub, repeat(fexp), repeat(ppi_matrix), repeat(cellpear), repeat(k), repeat(lgene_exp), cells, repeat(maxSR))

    SRs = []
    for i in ans:
        SRs.append(i)
    with open(result_file,'w') as f:
        SRst = [str(tp) for tp in SRs]
        f.writelines('\n'.join(SRst))
