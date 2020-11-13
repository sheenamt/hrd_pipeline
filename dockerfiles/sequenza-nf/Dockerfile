FROM "sequenza/sequenza:3.0.0"

WORKDIR /
USER root
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
	&& (apt-get install -y procps) \
	&& apt-get clean && apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*

# update sequenza-utils from Pypi
RUN wget https://files.pythonhosted.org/packages/86/b7/5343a87e787affca67fa02097259e90d15eb9c0f1ffab3972b8193958389/sequenza-utils-3.0.0.tar.gz -O sequenza-utils.tar.gz \
    && tar -xvpf sequenza-utils.tar.gz \
    && cd sequenza-utils* \
    && python setup.py install --install-scripts=/usr/bin \
    && cd ../ && rm -rf *sequenza-utils*