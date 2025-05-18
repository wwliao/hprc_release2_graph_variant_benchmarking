FROM mambaorg/micromamba:1.5.10-noble

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER calculate_metrics.py ./
RUN chmod +x /app/calculate_metrics.py

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        cyvcf2=0.31.1 procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y

USER root
ENV PATH="/app:$MAMBA_ROOT_PREFIX/bin:$PATH"
