# Use the official Ubuntu base image
#FROM ubuntu:latest
# Set environment variables
# ENV DEBIAN_FRONTEND noninteractive

# Install system dependencies
#RUN apt-get update && apt-get install -y \
#    python3.11 \
#    python3.11-dev \
#    python3-pip \
#    && rm -rf /var/lib/apt/lists/*

FROM python:3.9-slim

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install PoreAnalyser

# Set the working directory in the container
WORKDIR /app

# Copy the requirements.txt file into the container
# COPY requirements.txt .
#COPY entrypoint.sh /app/entrypoint.sh
COPY app.py /app/app.py
COPY PoreAnalyser /app/PoreAnalyser
COPY visualise_pathway_hole.tcl /app/visualise_pathway_hole.tcl
COPY chimera_pore.py /app/chimera_pore.py
COPY pymol_pore_visu.py /app/pymol_pore_visu.py
COPY README.md /app/README.md
COPY chimeraX_pore.py /app/chimeraX_pore.py

# Install Python dependencies
#RUN python3.11 -m pip install PoreAnalyser #--no-cache-dir -r requirements.txt

# The EXPOSE instruction informs Docker that the container listens on the specified network ports at runtime. Your container needs to listen to Streamlit’s (default) port 8501:
EXPOSE 8501

# The HEALTHCHECK instruction tells Docker how to test a container to check that it is still working. Your container needs to listen to Streamlit’s (default) port 8501:
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

# Set the command to run your application
#CMD ["python3.11", "your_main_script.py"]

# Code file to execute when the docker container starts up (`entrypoint.sh`)
#ENTRYPOINT ["/entrypoint.sh"]
# An ENTRYPOINT allows you to configure a container that will run as an executable. Here, it also contains the entire streamlit run command for your app, so you don’t have to call it from the command line:
ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]

# https://docs.streamlit.io/knowledge-base/tutorials/deploy/docker

# Based on your server's network configuration, you could map to port 80/443 so that users can view your app using the server IP or hostname. For example: http://your-server-ip:80 or http://your-hostname:443.
