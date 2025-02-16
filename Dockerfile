FROM nfcore/base:2.1
# Dockerfile adapted from Felipe https://github.com/fmalmeida/ngs-preprocess/blob/master/Dockerfile
LABEL authors="Henrietta Holze" \
      description="Docker image containing all software requirements for the BARtab nextflow pipeline"

# Install the conda environment
RUN conda install -y -c conda-forge mamba
COPY environment.yaml /
RUN mamba env create --quiet -f /environment.yaml && mamba clean -a
# remove fastp from environment.yaml and install afterwards 
# because mamba create never finishes when adding it to environment
RUN conda install --name bartab-1.4 -c bioconda fastp=0.23.4

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/bartab-1.4/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name bartab-1.4 > bartab-1.4.yaml
