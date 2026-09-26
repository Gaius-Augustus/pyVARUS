# Manual installation of pyVARUS

The [container](../README.md#container-recommended) and the
[conda environment](../README.md#conda) cover most users. This page is for
installing without either, e.g. with your distro's package manager.

pyVARUS has two kinds of dependencies: Python packages (installed by `pip`)
and external command-line tools (installed by you).

## 1. External command-line tools

These are *not* installed by `pip` and must be on `PATH` before you run pyVARUS.

| Tool | Used by | Required? |
|---|---|---|
| `hisat2`, `hisat2-build` | `varus index`, `varus run` (short reads) | required unless using `--longreads` |
| `minimap2` | `varus run` (Logan pre-screen), `varus logan`, `varus index --longreads`, `varus run --longreads` | required unless using `--no-logan` |
| `samtools` (>= 1.17) | `varus run` (sort, merge, index) | required |
| `fastq-dump` ([sra-toolkit](https://github.com/ncbi/sra-tools)) | `varus run` (downloads from SRA) | required |
| `zstd` | Logan contig decompression, only if the `zstandard` Python package is missing | optional |

References for all tools are listed under [Citation](../README.md#citation).

Install via conda into an existing environment:

```sh
conda install -c conda-forge -c bioconda hisat2 minimap2 "samtools>=1.17" sra-tools zstd
```

Or via your distro package manager (Ubuntu example; check that the packaged
samtools is at least 1.17):

```sh
sudo apt install hisat2 minimap2 samtools sra-toolkit zstd
```

Then disable the NCBI cache as described in the
[README](../README.md#disable-the-ncbi-cache).

## 2. Python package

```sh
git clone https://github.com/Gaius-Augustus/pyVARUS.git
cd pyVARUS
pip install -e ".[align]"   # add ',dev' for the test suite
```

The `[align]` extra pulls in `pysam` (needed by `varus run` for intron
extraction). `pip` builds it from source against `htslib`, which only
compiles on Linux/macOS -- on Windows you can still install plain
`pip install -e .` to use the `runlist` and `index` subcommands.
