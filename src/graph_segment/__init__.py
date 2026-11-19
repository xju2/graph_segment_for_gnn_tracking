"""Graph segmentation for GNN tracking: one entry point, `get_tracks`, over interchangeable backends.

    tracks = get_tracks(num_nodes, edge_index, scores, "cc_bestpred_walk", th_min=0.3, th_add=0.8)

Every backend takes the same inputs -- edge_index (2, E) pointing outward (from the hit
nearer the interaction point, so the graph is acyclic) and scores (E,) -- and returns
{method: [[node, ...], ...]}, so switching backend is changing one string.
To add a backend: write `fn(num_nodes, edge_index, scores, **cuts) -> dict` and put it in BACKENDS.
"""

from __future__ import annotations

import inspect

from . import _core
from ._core import cc_and_bestpred_walk

__all__ = ["BACKENDS", "cc_and_bestpred_walk", "get_tracks"]


def _cc_walk(num_nodes, edge_index, scores, walk, best_predecessor, cc_cut, th_min=0.3,
             th_add=0.8, min_track_score=0.0, min_track_hits=0):
    cc, w = _core.cc_and_walk(num_nodes, edge_index, scores, cc_cut, th_min, th_add, walk,
                              best_predecessor, min_track_score, min_track_hits)
    return {"cc": cc, "walk": w} if walk else {"cc": cc}


def cc(num_nodes, edge_index, scores, cc_cut=0.01, min_track_score=0.0, min_track_hits=0):
    """Connected components after cc_cut; keep the ones that are simple paths."""
    return _cc_walk(num_nodes, edge_index, scores, False, False, cc_cut, 0.0, 0.0,
                    min_track_score, min_track_hits)


def cc_walk(num_nodes, edge_index, scores, cc_cut=0.01, th_min=0.3, th_add=0.8,
            min_track_score=0.0, min_track_hits=0):
    """cc, then walk through the branching components: from each start hit follow the best
    edge above th_min and branch on every edge above th_add, keep the longest path."""
    return _cc_walk(num_nodes, edge_index, scores, True, False, cc_cut, th_min, th_add,
                    min_track_score, min_track_hits)


def cc_bestpred_walk(num_nodes, edge_index, scores, cc_cut=0.01, th_min=0.3, th_add=0.8,
                     min_track_score=0.0, min_track_hits=0):
    """cc_walk, stepping u -> v only if u -> v is v's best edge from an unused hit."""
    return _cc_walk(num_nodes, edge_index, scores, True, True, cc_cut, th_min, th_add,
                    min_track_score, min_track_hits)


def cc_junction_removal(num_nodes, edge_index, scores, cc_cut=0.7, junction_cut=0.5,
                        min_chain_length=3, min_track_score=0.0, min_track_hits=0):
    """Simple chains after cc_cut; on the other hits keep each edge above junction_cut that
    is the best in-edge of its target and the best out-edge of its source."""
    cc_, jr = _core.cc_and_junction_removal(num_nodes, edge_index, scores, cc_cut, junction_cut,
                                            min_chain_length, min_track_score, min_track_hits)
    return {"cc": cc_, "junction": jr}


BACKENDS = {
    "cc": cc,
    "cc_walk": cc_walk,
    "cc_bestpred_walk": cc_bestpred_walk,
    "cc_junction_removal": cc_junction_removal,
}

_KNOWN = {p for f in BACKENDS.values() for p in inspect.signature(f).parameters}


def get_tracks(num_nodes, edge_index, scores, backend="cc_bestpred_walk", **cuts):
    """
    Segment a scored graph into tracks with the named backend (see BACKENDS).

    cuts may hold the cuts of several backends; each backend takes the ones it knows,
    and a name no backend knows (a typo) raises. Every backend takes min_track_score and
    min_track_hits: drop tracks whose mean internal edge score is below min_track_score,
    or with fewer than min_track_hits hits.
    Returns {method: [[node, ...], ...]}, node indices into the graph.
    """
    if backend not in BACKENDS:
        msg = f"unknown backend {backend!r}, choose from {list(BACKENDS)}"
        raise ValueError(msg)
    unknown = set(cuts) - _KNOWN
    if unknown:
        msg = f"unknown segmentation cuts {sorted(unknown)}"
        raise TypeError(msg)
    fn = BACKENDS[backend]
    accepted = inspect.signature(fn).parameters
    return fn(num_nodes, edge_index, scores, **{k: v for k, v in cuts.items() if k in accepted})
