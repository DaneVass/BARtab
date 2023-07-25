FROM nfcore/base
# Dockerfile adapted from Felipe https://github.com/fmalmeida/ngs-preprocess/blob/master/Dockerfile
LABEL authors="Henrietta Holze" \
      description="Docker image containing all software requirements for the BARtab nextflow pipeline"

# Install the conda environment
RUN conda install -y -c conda-forge mamba
COPY environment.yml /
RUN mamba env create --quiet -f /environment.yml && mamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/bartab-1.3/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name bartab-1.3 > bartab-1.3.yml
