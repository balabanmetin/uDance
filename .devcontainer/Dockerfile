FROM continuumio/miniconda3:latest
RUN apt-get update -qq && apt-get upgrade -qq -y && \
    apt-get install -qq -y default-jdk && \
    git clone https://github.com/balabanmetin/uDance.git && \
    cd uDance && \
    bash install.sh && \
    conda activate udance && \
    echo "conda activate udance" >> ~/.bashrc
ENTRYPOINT ["bash", "-l", "-c"]
