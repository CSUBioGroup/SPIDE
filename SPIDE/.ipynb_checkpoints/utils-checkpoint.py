import numpy as np
import networkx as nx
import pandas as pd

def extract_ppi(ppifile):
    ppi_df = pd.read_csv(ppifile)
    gene_id = {}
    id_gene = {}
    for idx, row in ppi_df.iterrows():
        g1, g2 = int(row['gene1']), int(row['gene2'])
        if g1 in gene_id.keys():
            continue
        else:
            newid = len(gene_id)
            gene_id[g1] = newid
            id_gene[newid] = g1

        if g2 in gene_id.keys():
            continue
        else:
            newid = len(gene_id)
            gene_id[g2] = newid
            id_gene[newid] = g2
    gene_num = len(gene_id)
    ppi_matrix = np.zeros((gene_num, gene_num))
    for idx, row in ppi_df.iterrows():
        g1, g2 = int(row['gene1']), int(row['gene2'])
        g1id = gene_id[g1]
        g2id = gene_id[g2]
        ppi_matrix[g1id][g2id] = 1
        ppi_matrix[g2id][g1id] = 1
    
    return ppi_matrix, gene_id, id_gene

def DoIntegPPI(gene_exp, ppi_matrix, gene_id, id_gene):
    inter_genes = set(gene_exp.keys())&set(gene_id.keys())
    
    gene_num = len(gene_id)
    inter_matrix = []
    for idx in range(gene_num):
        if id_gene[idx] in inter_genes:
            continue
        else:
            ppi_matrix[:, idx] = 0
            ppi_matrix[idx, :] = 0
    ppi_matrix = np.maximum(ppi_matrix,ppi_matrix.T)
    g = nx.from_numpy_array(ppi_matrix)
    largest_components=max(nx.connected_components(g),key=len)
    subg = g.subgraph(largest_components)
    
    new_old = {}
    old_new = {}
    new_gene_id = {}
    new_id_gene = {}
    for newidx, oldidx in enumerate(subg.nodes()):
        new_old[newidx] = oldidx
        old_new[oldidx] = newidx
        gname = id_gene[oldidx]
        new_gene_id[gname] = newidx
        new_id_gene[newidx] = gname
    
    new_ppi = np.zeros((len(largest_components), len(largest_components)))
    for idx in range(new_ppi.shape[0]):
        new_ppi[idx][idx] = 1.0
    for edge_s, edge_t in subg.edges():
        s = old_new[edge_s]
        t = old_new[edge_t]
        new_ppi[s][t] = 1.0
        new_ppi[t][s] = 1.0
    
    new_gene_exp = []
    for idx in range(new_ppi.shape[0]):
        gname = new_id_gene[idx]
        new_gene_exp.append(gene_exp[gname])
    new_gene_exp = np.array(new_gene_exp)
    
    return new_gene_exp, new_ppi, new_gene_id, new_id_gene


def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

