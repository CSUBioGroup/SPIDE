Directly Test SPIDE (without install)
```
cd SPIDE
python pipeline.py
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