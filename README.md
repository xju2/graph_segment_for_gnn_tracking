# Graph Segment Algorithms for Graph Neural Network-based Tracking

### Installation
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

### Usage
Download the data `python3 download_data.py`. And then `./build/bin/walk_through data/debug_graph.dot`.

```text
Input Graph: 271663 vertices, 615710 edges.
From CC&&Walk: Number of tracks found by Walkthrough: 1260
From CC&&Walk:: Total 4604 tracks
Time taken: 186.526 ms
From ACORN: Number of tracks found by CC: 2949
From ACORN: Number of tracks found by Walkthrough: 1299
From ACORN: Total 4248 tracks.
```
