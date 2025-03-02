# Graph Segment Algorithms for Graph Neural Network-based Tracking

### Installation
```bash
pip install .
```

### Test
Download the data `python3 download_data.py`. And then `./build/bin/walk_through data/debug_graph.dot`.


### Developer Guide
The following commands are for Perlmutter. 
For other systems, please adjust the commands accordingly.
```bash
podman-hpc pull docker.io/docexoty/mltools:20250227
```
Or you can build the image from the [Dockerfile](Dockerfile) in this repository.

Then run the following command to start the container,
and build the code:
```bash
podman-hpc run -it --rm --gpu -v $PWD:$PWD -w $PWD docexoty/mltools:20250227 bash
```
```bash
cmake -B build -S . -Dpybind11_DIR=/usr/local/lib/python3.10/dist-packages/pybind11/share/cmake/pybind11
cmake --build build
```

### Packing and Uploading
```bash
python3 -m build --sdist
python3 -m build --wheel
```
#### Note
The code was improved by `OpenAI o1Pro` model. The execution time was reduced from 403 ms to 190 ms.
