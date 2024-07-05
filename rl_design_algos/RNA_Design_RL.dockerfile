# Use the official Debian image as the base
FROM debian:latest

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    git \
    gcc \
    g++ \
    && apt-get clean

# Download and install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Update the PATH environment variable
ENV PATH /opt/conda/bin:$PATH

# Set the working directory
WORKDIR /app

# Copy the environment file
COPY environment.yml /app/

# Create the Conda environment
RUN conda env create -f environment.yml
# RUN echo "source activate rna_design_rl" > ~/.bashrc
RUN /bin/bash -c "source activate rna_design_rl"
ENV PATH /opt/conda/envs/myenv/bin:$PATH
ENV CONDA_DEFAULT_ENV rna_design_rl

# Activate the environment
# SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Copy the rest of the application code
COPY . /app/

# Specify the default command to run the main.py script
CMD ["conda", "run", "-n", "rna_design_rl", "python", "main.py"]
