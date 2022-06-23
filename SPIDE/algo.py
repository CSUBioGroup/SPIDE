import numpy as np
import pandas as pd
from tqdm import tqdm

def fillfunc2(gene_exp, cellpear, k):
    fill_gene_exp = gene_exp.copy()
    cell_num = cellpear.shape[0]
    
    for cell in tqdm(range(cell_num)):
        sortidx = np.argsort(cellpear[cell])
        knn_exp = []
        for i in range(k):
            cellidx = sortidx[cell_num-i-1]
            cell_exp = gene_exp[:, cellidx].tolist()
            knn_exp.append(cell_exp)
        knn_exp = np.array(knn_exp)
        knn_mean = np.mean(knn_exp, axis=0)
        back_mean = np.mean(gene_exp.T, axis=0)
        
        count = 0 
        for idx, val in enumerate(fill_gene_exp[:,cell]):
            if val==0 and knn_mean[idx]>back_mean[idx]:
                count += 1
                fill_gene_exp[idx,cell] = knn_mean[idx]
    return fill_gene_exp


def applycomps(p_matrix):
    gene_num = p_matrix.shape[0]
    S_v = [0]*gene_num
    for idx in range(gene_num):
        for val in p_matrix[:, idx]:
            if val>0:
                S_v[idx] -= val*np.log(val)
    S_v = np.array(S_v)
    return S_v

def newSCENT2ps_abs_fill(gene_exp, ppi_matrix, cellpear, k, lgene_exp, local=True, maxSR=0):
    SR_vals = []
    cell_num = cellpear.shape[0]
    gene_num = gene_exp.shape[0]
    for cell in tqdm(range(cell_num)):
        sortidx = np.argsort(cellpear[cell])
        knn_exp = []
        for i in range(k):
            cellidx = sortidx[cell_num-i-1]
            cell_exp = lgene_exp[:, cellidx].tolist()
            knn_exp.append(cell_exp)
        knn_exp = np.array(knn_exp) # cell-gene
        knn_exp = knn_exp.T # gene-cell
        genepear = np.corrcoef(knn_exp)

        genepear[np.isnan(genepear)] = 0.001
        genepear[np.isinf(genepear)] = 0.001
        weighted_ppi_matrix = np.abs(ppi_matrix*genepear)
        # weighted_ppi_matrix[np.isnan(weighted_ppi_matrix)] = 0.001
        row, col = np.diag_indices_from(weighted_ppi_matrix)
        weighted_ppi_matrix[row,col] = 0

        exp_v = gene_exp[:, cell]
        sumexp_v = np.dot(weighted_ppi_matrix, exp_v)
        invP_v = exp_v*sumexp_v
        nf = np.sum(invP_v)
        invP_v = invP_v/nf

        p_m =(weighted_ppi_matrix*exp_v).T/sumexp_v
        S_v = applycomps(p_m)
        
        SR = np.sum(invP_v*S_v)
        
        if maxSR>0:
            SR = SR/maxSR
        SR_vals.append(SR)
        print(cell, SR)
    return SR_vals

def newSCENT2ps_abs_fill_sub(gene_exp, ppi_matrix, cellpear, k, lgene_exp, cell, local=True, maxSR=0):
    SR_vals = []
    cell_num = cellpear.shape[0]
    gene_num = gene_exp.shape[0]
    
    sortidx = np.argsort(cellpear[cell])
    knn_exp = []
    for i in range(k):
        cellidx = sortidx[cell_num-i-1]
        cell_exp = lgene_exp[:, cellidx].tolist()
        knn_exp.append(cell_exp)
    knn_exp = np.array(knn_exp) # cell-gene
    knn_exp = knn_exp.T # gene-cell
    genepear = np.corrcoef(knn_exp)

    genepear[np.isnan(genepear)] = 0.001
    genepear[np.isinf(genepear)] = 0.001
    weighted_ppi_matrix = np.abs(ppi_matrix*genepear)
    # weighted_ppi_matrix[np.isnan(weighted_ppi_matrix)] = 0.001
    row, col = np.diag_indices_from(weighted_ppi_matrix)
    weighted_ppi_matrix[row,col] = 0

    exp_v = gene_exp[:, cell]
    sumexp_v = np.dot(weighted_ppi_matrix, exp_v)
    invP_v = exp_v*sumexp_v
    nf = np.sum(invP_v)
    invP_v = invP_v/nf

    p_m =(weighted_ppi_matrix*exp_v).T/sumexp_v
    S_v = applycomps(p_m)

    SR = np.sum(invP_v*S_v)
    if maxSR>0:
        SR = SR/maxSR
    return SR