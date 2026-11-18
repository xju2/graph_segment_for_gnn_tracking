# Graph Segment Algorithms for Graph Neural Network-based Tracking

### Installation
```bash
pip install graph_segment
```

### CC + best-predecessor walk
`cc_and_bestpred_walk` (`src/cpp_core/bestpred_walk.cxx`) works on plain edge arrays, no Boost.
The walk moves from hit u to hit v only if u -> v is v's best incoming edge among the hits not
yet used, so a branch that joins another track's chain cannot take its hits. Optional
post-filters drop tracks by mean edge score and number of hits.
```python
import graph_segment
# edge_index: (2, E) int array, pointing outward (from the hit nearer the IP), scores: (E,) floats
cc_tracks, walk_tracks = graph_segment.cc_and_bestpred_walk(
    num_nodes, edge_index, scores, cc_cut=0.01, th_min=0.3, th_add=0.8,
    min_track_score=0.5, min_track_hits=7)   # lists of node-index lists
```
It gives the same tracks as acorn's `CCandWalk` with `bestpred_walk: True` (checked on 993 muon
collider events; `tests/test_bestpred_walk.py` compares the two on random graphs when acorn is
installed). Run the tests with `pytest tests/`.

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
