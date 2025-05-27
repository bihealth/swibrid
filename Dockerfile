# Use the miniconda base image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Create a non-root user and group
RUN groupadd -r swibridgroup && useradd -r -g swibridgroup swibriduser

# Set umask to ensure correct permissions
RUN echo "umask 002" >> /etc/profile

# Create directories and set ownership
RUN mkdir -p /home/swibriduser/.cache && chown -R swibriduser:swibridgroup /home/swibriduser

# Copy the environment file to the Docker image
COPY swibrid_env.yaml .

# Create the conda environment
RUN conda env create -f swibrid_env.yaml

# Setup SSH keys 
RUN mkdir -p /root/.ssh
COPY id_rsa /root/.ssh/id_rsa
RUN chmod 600 /root/.ssh/id_rsa

# Add GitHub to known hosts to prevent SSH from asking for confirmation
RUN ssh-keyscan github.com >> /root/.ssh/known_hosts

# Clone the swibrid repository 
RUN git clone git@github.com:bihealth/swibrid.git

# Install swibrid
WORKDIR /app/swibrid
RUN conda run -n swibrid_env python setup.py install

# Set environment variables
ENV PATH=/opt/conda/envs/swibrid_env/bin:$PATH
ENV XDG_CACHE_HOME=/home/swibriduser/.cache

# Set the default command to activate the environment
CMD ["conda", "run", "-n", "swibrid_env", "/bin/bash"]

# Switch to the non-root user
USER swibriduser

# re-set the working dir
WORKDIR /home/swibriduser/ 
ENTRYPOINT ["swibrid"]
