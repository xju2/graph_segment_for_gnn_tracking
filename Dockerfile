# docexoty/mltools:20250227
# BASE image: https://hub.docker.com/r/nvidia/cuda
FROM nvidia/cuda:12.4.1-cudnn-devel-ubuntu22.04

RUN apt-get update -y && DEBIAN_FRONTEND="noninteractive" TZ="America/Los_Angeles" apt-get -y install \
	dpkg-dev g++ gcc python3-dev python3-pip \
    libboost-all-dev vim wget rsync \
    libtbb2 libtbb-dev cmake \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip && pip install pybind11
RUN pip install cugraph-cu12 nx-cugraph-cu12 cugraph-pyg-cu12 --extra-index-url=https://pypi.nvidia.com
