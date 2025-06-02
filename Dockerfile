FROM ubuntu:20.04 AS base

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHON_VERSION=3.8.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential curl wget git vim \
    libssl-dev zlib1g-dev libbz2-dev libreadline-dev \
    libsqlite3-dev llvm libncurses5-dev libncursesw5-dev \
    xz-utils tk-dev libffi-dev liblzma-dev \
    python3-pip python3.8 python3.8-venv python3-dev \
    r-base r-base-dev ghostscript \
    slurm-wlm munge libmunge-dev rng-tools \
    && apt-get clean

FROM base AS python-deps

# Install Python packages in a venv
RUN python3.8 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --upgrade pip && \
    pip install \
    pandas==1.3.0 \
    numpy==1.21.0 \
    networkx==2.2 \
    scikit-learn==0.24.2 \
    matplotlib==3.5.1 \
    seaborn==0.11.2 \
	scikit-learn-extra==0.2.0

FROM python-deps AS r-deps

# Install remotes for version-specific package install
RUN Rscript -e "install.packages(c('remotes'), repos='https://cloud.r-project.org')"

# Install specific versions of R packages
RUN Rscript -e "remotes::install_version('FMStable', version = '0.1-4', repos = 'https://cloud.r-project.org')" && \
    Rscript -e "remotes::install_version('harmonicmeanp', version = '3.0.1', repos = 'https://cloud.r-project.org')"

FROM r-deps AS slurm-setup

# Set up Munge
RUN mkdir -p /etc/munge /var/run/munge /var/log/munge && \
    chown -R munge:munge /etc/munge /var/run/munge /var/log/munge && \
    chmod 0700 /etc/munge /var/run/munge /var/log/munge && \
    /usr/sbin/create-munge-key && \
    chown munge:munge /etc/munge/munge.key

# SLURM minimal config
COPY slurm.conf /etc/slurm-llnl/slurm.conf
RUN chmod 644 /etc/slurm-llnl/slurm.conf
RUN mkdir -p /var/spool/slurmctld /var/spool/slurmd /var/log/slurm && \
    chown -R slurm:slurm /var/spool/slurmctld /var/spool/slurmd /var/log/slurm

FROM slurm-setup AS final

# Add code base
COPY . /doublethink-mcmc

# Set working directory
RUN mkdir -p /workspace
WORKDIR /workspace

# Start munge, slurmctld, and slurmd in the background
COPY start_slurm.sh /usr/local/bin/start_slurm.sh
RUN chmod +x /usr/local/bin/start_slurm.sh

# Execute start-up script
CMD ["/usr/local/bin/start_slurm.sh"]
