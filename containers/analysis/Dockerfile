# Python analysis image for all scripts in bin/.
# Built and published to GHCR by .github/workflows/container.yml; Nextflow pins it by tag.
FROM mambaorg/micromamba:2.10.0-ubuntu24.04
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba list -n base > /opt/conda/conda-list.txt && \
    micromamba clean -a -y
# Apptainer ignores the micromamba entrypoint, so put the env on PATH explicitly.
ENV PATH=/opt/conda/bin:$PATH
