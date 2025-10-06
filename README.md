# Graph Segment Algorithms for Graph Neural Network-based Tracking

### Installation
```bash
pip install graph_segment
```

### Test
Download the data `python3 download_data.py`. And then `./build/bin/walk_through data/debug_graph.dot`.


### Developer Guide
Install a Python virtual environment and activate it:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install cibuildwheel build twine
```

Then you can run
```bash
pipx run cibuildwheel --platform linux
```

### Packing and Uploading
```bash
python3 -m build --sdist
twine upload dist/* wheelhouse/*
```
#### Note
The code was improved by `OpenAI o1Pro` model. The execution time was reduced from 403 ms to 190 ms.
