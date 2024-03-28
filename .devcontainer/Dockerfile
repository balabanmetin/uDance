FROM continuumio/miniconda3:latest
RUN git clone https://github.com/balabanmetin/uDance.git && \
    cd uDance && \
    bash install.sh && \
    echo "conda activate udance" >> ~/.bashrc
ENTRYPOINT ["bash", "-l", "-c"]
