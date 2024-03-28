FROM continuumio/miniconda3:latest
RUN conda install -y openjdk && \
    git clone https://github.com/balabanmetin/uDance.git && \
    cd uDance && \
    bash install.sh && \
    conda activate udance && \
    echo "conda activate udance" >> ~/.bashrc
ENTRYPOINT ["bash", "-l", "-c"]
