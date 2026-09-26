# pyVARUS: online sampling of complementary RNA-seq reads from NCBI SRA
# (and Logan contigs) for genome annotation.
#
# Includes: pyVARUS [align], HISAT2, minimap2, samtools, sra-tools, zstd
#
# The image is built from the checked-out source tree, not from a git clone,
# so it contains exactly the commit you build from. Build context is the
# repository root:
#
#   docker/build_push.sh                      # build + push, tag from pyproject
#   docker build -f docker/Dockerfile -t gaiusaugustus/pyvarus:<tag> .
#
# Convert to Singularity:
#   singularity pull pyvarus.sif docker://gaiusaugustus/pyvarus:latest
#
# Usage:
#   varus runlist "Genus species" --outdir /out --email user@host
#   varus index genome.fa --outdir /out --threads 8
#   varus run "Genus species" genome.fa --runlist /out/runlist.tsv \
#       --index /out/genome_index --outdir /out --threads 8

FROM mambaorg/micromamba:1.5.8

USER root

# Bioconda tools. samtools>=1.17 for CSI index support (-c flag).
# zstd is the fallback decompressor for Logan contigs if the zstandard
# Python package fails to import.
RUN micromamba install -y -n base \
        -c conda-forge -c bioconda \
        python=3.11 \
        hisat2 \
        minimap2 \
        "samtools>=1.17" \
        sra-tools \
        zstd \
        pip && \
    micromamba clean -afy

ENV PATH="/opt/conda/bin:$PATH"

# pyVARUS with the [align] extra (pysam for intron extraction from BAM).
COPY pyproject.toml README.md LICENSE /tmp/pyvarus/
COPY varus /tmp/pyvarus/varus
RUN pip install --no-cache-dir "/tmp/pyvarus[align]" && \
    rm -rf /tmp/pyvarus

ARG VARUS_GIT_COMMIT=unknown
LABEL org.opencontainers.image.source="https://github.com/Gaius-Augustus/pyVARUS" \
      org.opencontainers.image.revision="${VARUS_GIT_COMMIT}"

# Smoke test: entry point and every external binary pyVARUS calls.
RUN varus --help > /dev/null && \
    python -c "import pysam, zstandard, pyfaidx, Bio, numpy" && \
    for b in hisat2 hisat2-build minimap2 samtools fastq-dump zstd; do \
        command -v "$b" > /dev/null || { echo "missing: $b"; exit 1; }; \
    done

ENTRYPOINT ["varus"]
CMD ["--help"]
