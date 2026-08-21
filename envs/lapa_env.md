# lapa_env — LAPA 0.0.5 conda env recipe

LAPA 0.0.5's pip dependencies (ncls / sorted-nearest / bamread) are old Cython
extensions that FAIL to compile under Cython 3, and `pip install lapa` also silently
falls back to a broken system-python install if the env isn't active. The working
recipe pins `cython<3`, gets `ncls` from bioconda, disables build isolation, and pins
pyranges to the 0.0.x line LAPA was written against.

```bash
ENV=/sc/arion/projects/als-omics/conda/envs/lapa_env
export CONDA_PKGS_DIRS=<a-writable-scratch>/.conda_pkgs   # shared pkgs cache isn't writable

mamba create -y -p $ENV -c conda-forge -c bioconda \
  python=3.9 "numpy<=1.23" pandas scipy pybigwig pysam samtools matplotlib \
  click tqdm "cython<3" pip setuptools wheel ncls
conda activate $ENV
pip install --no-build-isolation sorted-nearest pyranges bamread kipoiseq betabinomial
pip install --no-deps lapa
# LAPA 0.0.5 targets the pyranges 0.0.x API; pip pulls 0.1.4 (a rewrite) -> pin it back:
pip install --no-build-isolation "pyranges<0.1" "sorted-nearest==0.0.33"
```

Verified working set: python 3.9, lapa 0.0.5, pyranges 0.0.133, sorted-nearest 0.0.33,
bamread 0.0.20, numpy<=1.23. Check: `python -c "import lapa; from lapa.result import
LapaResult; import pyranges as pr; pr.from_dict({'Chromosome':['chr1'],'Start':[1],
'End':[9]}).cluster()"` should run clean.

Consumed by `workflows/lapa_pipeline.smk` (config key `lapa_env: "lapa_env"`).
