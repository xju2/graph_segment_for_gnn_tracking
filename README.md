# Graph Segment Algorithms for Graph Neural Network-based Tracking

### Installation
```bash
pip install graph_segment
```

### Usage: `get_tracks`
After a GNN has scored the edges of a hit graph, graph segmentation cuts the graph into
tracks. All algorithms ("backends") sit behind one function, so switching algorithm means
changing one string:
```python
import graph_segment

# edge_index: (2, E) int array, each edge pointing outward (from the hit nearer the
#             interaction point), so the graph is acyclic; scores: (E,) edge scores
tracks = graph_segment.get_tracks(num_nodes, edge_index, scores, "cc_bestpred_walk",
                                  cc_cut=0.01, th_min=0.3, th_add=0.8,
                                  min_track_score=0.5, min_track_hits=7)
# {"cc": [[node, ...], ...], "walk": [...]}   node indices into the graph
print(list(graph_segment.BACKENDS))
```

| backend | what it does |
|---|---|
| `cc` | drop edges with score <= `cc_cut`; the connected components that are simple paths (no hit with two incoming or two outgoing edges, >= 3 hits) are tracks |
| `cc_walk` | `cc`, then walk through the branching components. From each unused hit in topological order, follow the best edge above `th_min` and branch on every edge above `th_add`, then keep the longest path |
| `cc_bestpred_walk` | `cc_walk` with the best-predecessor rule: step u -> v only if u -> v is v's best edge from an unused hit, so a branch that joins another track's chain cannot take its hits |
| `cc_junction_removal` | simple chains of >= `min_chain_length` hits after `cc_cut`. On the other hits, keep each edge above `junction_cut` that is the best in-edge of its target and the best out-edge of its source, and take their components as tracks |

- **Cuts:** `get_tracks` takes the cuts of every backend. Each backend uses only its own cuts, so one set of cuts works for all of them, and a misspelled cut name raises an error.
- **Filters:** every backend takes `min_track_score` (drop tracks whose mean internal edge score is lower) and `min_track_hits`.
- **Track length:** no backend returns a track with fewer than two hits.
- **Implementation:** the backends are C++ on plain edge arrays (`src/cpp_core/bestpred_walk.cxx`), with no Boost and no networkx.
- **Adding a backend:** write `fn(num_nodes, edge_index, scores, **cuts) -> {method: tracks}` and add it to `BACKENDS` in `src/graph_segment/__init__.py`.

`cc_and_bestpred_walk(...)` is still available and returns `(cc_tracks, walk_tracks)`.

**Validation:** every backend gives the same tracks as the Python implementation in
acorn (`acorn.stages.track_building.segmentation`). The one exception is that acorn's
plain walk also returns dead-end hits as 1-hit "tracks". `tests/test_get_tracks.py` makes
this comparison on random graphs when acorn is installed. Run the tests with `pytest tests/`.

### Benchmark: muon collider with BIB
Setup:
- **Events:** 993 validation events, each with one muon plus ~120k beam-induced-background
  hits (~72k GNN edges per event). The GNN is acorn's v3 muon-collider model (epoch 154).
- **Script:** `benchmarks/muon_collider.py`. It needs acorn to read the graphs and match tracks.
- **Matching:** acorn's matching with 50% of hits and tracks of >= 5 hits.
- **Metrics:** *eff* uses ATLAS matching and *eff 2-way* two-way matching. *fake* is the
  fraction of tracks matched to no particle. *dup* is the number of extra matched tracks per
  reconstructed muon.
- **Filtered:** `min_track_score=0.5`, `min_track_hits=7`.
- **Time:** `get_tracks` wall time per event with one thread on a Perlmutter CPU node.
  Loading the graph and matching are excluded.

| backend | cuts | eff | fake | dup | eff 2-way | eff / fake / dup (filtered) | time / event |
|---|---|---|---|---|---|---|---|
| `cc` | 0.01 | 0.000 | 1.000 | - | 0.000 | 0.000 / - / - | 9.6 ms |
| `cc` | 0.9 | 0.189 | 0.041 | 0.000 | 0.178 | 0.183 / 0.005 / 0.000 | 10.2 ms |
| `cc_walk` | 0.01 / 0.3 / 0.8 | 0.989 | 0.306 | 0.017 | 0.983 | 0.985 / 0.027 / 0.005 | 10.7 ms |
| **`cc_bestpred_walk`** | 0.01 / 0.3 / 0.8 | **0.993** | 0.303 | 0.012 | **0.988** | **0.990 / 0.019 / 0.003** | **10.2 ms** |
| `cc_junction_removal` | 0.7 / 0.5 / 3 hits | 0.991 | 0.034 | 0.138 | 0.928 | 0.979 / 0.003 / 0.052 | 16.4 ms |

The physics numbers are identical to acorn's Python backends on the same events. The speed
is not: acorn takes 1.15 s per event for `cc` and `cc_walk`, which use networkx, and about
22 ms for its numpy `cc_bestpred_walk` and `cc_junction_removal`.

What to take from it:
- **Use `cc_bestpred_walk` with the filter.** It has the best efficiency with a 2% fake rate.
- **Plain CC fails with this background.** At a loose cut every muon component touches
  background hits and branches. At a tight cut the muon breaks into pieces.
- **Without the filter every walk makes ~30% fakes.** These are mostly short 5–6 hit
  background segments. The length and score filter removes them for < 0.5% efficiency.
- **Junction removal splits the muon.** Its ATLAS efficiency and fake rate look good, but
  its dup rate is high and its 2-way efficiency is low. Compare backends on the 2-way
  efficiency and dup columns too.

Per-event data: https://portal.nersc.gov/cfs/m3443/xju/muon_collider/segmentation_benchmark/graph_segment/

### Test
`pytest tests/`


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
