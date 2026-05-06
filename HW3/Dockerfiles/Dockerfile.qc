FROM condaforge/mambaforge:latest
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge
RUN mamba install -y sra-tools fastqc trimmomatic entrez-direct && mamba clean -a -y
CMD ["/bin/bash"]
