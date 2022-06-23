# SPIDE

ppi file format (.csv):
```
gene1   gene2
10436   1510
7917    1510
1510    10436
...     ...
```

gene expression file format (.csv):
```
gene1, exp1, exp2, ..., expn
gene2, exp1, exp2, ..., expn
gene3, exp1, exp2, ..., expn
...
```

Directly Test SPIDE (without install)
```
python ./SPIDE/pipeline.py
```

Install SPIDE
```
pip install -e .
```

Import SPIDE
```
from SPIDE.spide import spide

spide(geneexp_file, k, ncore, ppi_file, result_file)
```

Details
```
spide(geneexp_file, k, ncore, ppi_file, result_file)

Parameters:
    geneexp_file: the path of single cell gene expression file
    k: the number of k-Nearest Neighbor
    ncore: the number of threads
    ppi_file: the path of underlying PPI file
    result_file: the path of result file (each row means the differential potency of corresponding cell)
```