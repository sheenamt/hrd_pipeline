FROM r-base:3.5.1
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       procps \
       libcurl4-openssl-dev libssl-dev \
       libxml2-dev \
       samtools \
       tabix \
       bwa \
       python python-dev python-setuptools python-pip \
       apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && wget https://files.pythonhosted.org/packages/86/b7/5343a87e787affca67fa02097259e90d15eb9c0f1ffab3972b8193958389/sequenza-utils-3.0.0.tar.gz -O sequenza-utils.tar.gz \
    && tar -xvpf sequenza-utils.tar.gz \
    && cd sequenza-utils* \
    && python setup.py install --install-scripts=/usr/bin \
    && cd ../ && rm -rf *sequenza-utils* \
    && mkdir /databases && chmod -R 7777 /databases \
    && mkdir /data && chmod -R 7777 /data \
    && wget https://bitbucket.org/ffavero/bio_pype/get/e723bf5.tar.gz -O bio_pype.tar.gz \
    && tar -xvpf bio_pype.tar.gz \
    && cd ffavero-bio_pype* && python setup.py test \
    && python setup.py install --install-scripts=/usr/bin \
    && cd ../ && rm -rf *bio_pype* \
    && pype repos install --force sequenza

VOLUME /databases /data

ADD exec/install_sequenza.R /usr/local/install_sequenza.R
RUN Rscript /usr/local/install_sequenza.R \
    && rm /usr/local/install_sequenza.R
ADD exec/run_sequenza.py /usr/bin/sequenza-pipeline
RUN chmod +x /usr/bin/sequenza-pipeline
RUN useradd -ms /bin/bash sequenza

USER sequenza
WORKDIR /home/sequenza

CMD ["/bin/bash"]


