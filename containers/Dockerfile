# All-in-one StrainQC-nf image — the registry twin of strainqc.def (same tool stack).
# Built + pushed to ghcr.io/digibio-ut/strainqc by .github/workflows/build-container.yml.
# Run the pipeline against it with:  -profile apptainer --strainqc_container docker://ghcr.io/digibio-ut/strainqc:<tag>
# (Apptainer/Docker both pull it; Apptainer converts to a SIF on first use.)
FROM mambaorg/micromamba:1.5.10

# gatk4 / bwa-mem2 pinned to match the validated hpc_modules + Apptainer runs; others float.
# procps-ng is required by Nextflow to collect task metrics (`ps`).
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        procps-ng \
        samtools \
        bwa-mem2=2.2.1 \
        gatk4=4.4.0.0 \
        breseq \
        kraken2 \
        mash \
        blast \
        cnvkit \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH" \
    LC_ALL=C \
    LANG=C

# Fail the build early if any tool is missing from PATH.
RUN for t in samtools bwa-mem2 gatk breseq kraken2 mash blastn cnvkit.py python3 ps; do \
        command -v "$t" >/dev/null || { echo "MISSING: $t"; exit 1; }; \
    done
