FROM r-base:3.5.1
# procps required for batchbot/nextflow
# -o APT::Immediate-Configure=0 fixes error <>

# Setup the base operating system and its requirements
RUN apt-get update \
    && apt-get install -y -o APT::Immediate-Configure=0 --no-install-recommends \
       libcurl4-openssl-dev libssl-dev \
       procps libgcc-8-dev gcc-8-base \
       libxml2-dev \
       samtools \
       tabix \
       bwa \
       python python-dev python-setuptools python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Setup the code for sequenza-utils
RUN wget https://files.pythonhosted.org/packages/86/b7/5343a87e787affca67fa02097259e90d15eb9c0f1ffab3972b8193958389/sequenza-utils-3.0.0.tar.gz -O sequenza-utils.tar.gz \
    && tar -xvpf sequenza-utils.tar.gz \
    && cd sequenza-utils* \
    && python setup.py install --install-scripts=/usr/bin \
    && cd ../ && rm -rf *sequenza-utils*

# Setup the code for sequenza 
RUN R -e "install.packages('sequenza',dependencies=TRUE, repos='http://cran.rstudio.com/')"

CMD ["/bin/bash"]
