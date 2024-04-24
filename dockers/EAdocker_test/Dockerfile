FROM ubuntu:20.04

LABEL version="0.5.0"
LABEL commit="setup antaRNA"

WORKDIR /app

# update apt and install essential requirements
RUN apt update && \
    apt install -y software-properties-common && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt update

# essential tools
RUN apt install -y wget tar build-essential vim

# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    /bin/bash miniconda.sh -b -p /opt/miniconda && \
    rm miniconda.sh
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc

ENV PATH="/opt/miniconda/bin:$PATH"
# create and active rnabench environment
RUN conda create -n rnabench python=3.9 -y

# --------------------------------------------
# |                antaRNA                   |
# --------------------------------------------
# ----------- package dependencies -----------
# ViennaRNA
WORKDIR /app
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_1_x/ViennaRNA-2.1.3.tar.gz && \
    tar -xzvf ViennaRNA-2.1.3.tar.gz && \
    rm ViennaRNA-2.1.3.tar.gz && \
    cd ViennaRNA-2.1.3 && \
    ./configure && \
    make && \
    make install
# ------------- python packages --------------
RUN conda install -n rnabench numpy scipy networkx -y
# ------------------ setup -------------------
COPY . /rnabench
RUN echo "conda activate rnabench && cd /rnabench/antarna && python setup.py install" >> ~/.bashrc
# RUN /bin/bash -c "source activate rnabench && cd /rnabench/antarna && python setup.py install"
# ----------- start -----------
RUN echo "conda activate rnabench" >> ~/.bashrc
WORKDIR /rnabench
CMD ["/bin/bash"]

