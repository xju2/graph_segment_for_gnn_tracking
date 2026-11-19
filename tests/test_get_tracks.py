from __future__ import annotations

import numpy as np
import pytest

from graph_segment import BACKENDS, get_tracks


def edges(*triples):
    e = np.array([t[:2] for t in triples]).T
    return e, np.array([t[2] for t in triples], dtype=np.float32)


# chain 0-1-2-3; hit 4 has a better edge into 2 than 1 has
STEAL = edges((0, 1, 0.9), (1, 2, 0.9), (2, 3, 0.9), (4, 2, 0.95))


def test_cc_only_keeps_simple_paths():
    e, s = edges((0, 1, 0.9), (1, 2, 0.9), (2, 3, 0.9), (5, 6, 0.9), (6, 7, 0.9), (5, 8, 0.9))
    assert get_tracks(9, e, s, "cc") == {"cc": [[0, 1, 2, 3]]}


def test_walk_with_and_without_best_predecessor():
    assert get_tracks(5, *STEAL, "cc_walk") == {"cc": [], "walk": [[0, 1, 2, 3]]}
    out = get_tracks(5, *STEAL, "cc_bestpred_walk")
    assert sorted(out["walk"]) == [[0, 1], [4, 2, 3]]


def test_junction_removal():
    # chain 0-1-2 above cc_cut; 3 -> 4 -> 5 and 6 -> 4 below it: 4's best in-edge is 3 -> 4
    e, s = edges((0, 1, 0.9), (1, 2, 0.9), (3, 4, 0.6), (4, 5, 0.6), (6, 4, 0.55))
    out = get_tracks(7, e, s, "cc_junction_removal", cc_cut=0.7, junction_cut=0.5,
                     min_chain_length=3)
    assert out == {"cc": [[0, 1, 2]], "junction": [[3, 4, 5]]}


def test_post_filters_for_every_backend():
    e, s = edges((0, 1, 0.9), (1, 2, 0.9), (2, 3, 0.9), (4, 5, 0.75), (5, 6, 0.75))
    for backend in BACKENDS:
        out = get_tracks(7, e, s, backend, cc_cut=0.7, min_chain_length=3,
                         min_track_score=0.8)
        assert [t for ts in out.values() for t in ts] == [[0, 1, 2, 3]], backend
        out = get_tracks(7, e, s, backend, cc_cut=0.7, min_chain_length=3, min_track_hits=4)
        assert [t for ts in out.values() for t in ts] == [[0, 1, 2, 3]], backend


def test_bad_arguments():
    with pytest.raises(ValueError, match="unknown backend"):
        get_tracks(5, *STEAL, "nope")
    with pytest.raises(TypeError, match="cc_cutt"):
        get_tracks(5, *STEAL, "cc", cc_cutt=0.1)


@pytest.mark.parametrize("backend", list(BACKENDS))
def test_matches_acorn_on_random_graphs(backend):
    """acorn's get_tracks (acorn.stages.track_building.segmentation) is the reference"""
    torch = pytest.importorskip("torch")
    seg = pytest.importorskip("acorn.stages.track_building.segmentation")
    from torch_geometric.data import Data

    cuts = {"cc_cut": 0.7, "junction_cut": 0.3, "min_chain_length": 3} \
        if backend == "cc_junction_removal" else {"cc_cut": 0.3}
    gen = torch.Generator().manual_seed(0)
    for _ in range(30):
        n, m = 300, 900
        ei = torch.randint(0, n, (2, m), generator=gen)
        ei = ei[:, ei[0] != ei[1]]
        g = Data(hit_id=torch.arange(n), hit_r=torch.rand(n, generator=gen),
                 hit_z=torch.rand(n, generator=gen), edge_index=ei,
                 edge_scores=torch.rand(ei.shape[1], generator=gen), num_nodes=n)
        ref = seg.get_tracks(g, backend, walk_min=0.4, walk_add=0.8, **cuts)  # points edges outward
        out = get_tracks(n, g.edge_index.numpy(), g.edge_scores.numpy(), backend,
                         th_min=0.4, th_add=0.8, **cuts)
        # acorn's plain walk also returns each dead-end start hit as a 1-hit "track"
        assert {frozenset(t) for ts in out.values() for t in ts} == \
            {frozenset(t) for ts in ref.values() for t in ts if len(t) > 1}


@pytest.mark.parametrize("backend", list(BACKENDS))
def test_no_single_hit_tracks(backend):
    rng = np.random.default_rng(0)
    e = rng.integers(0, 200, (2, 500))
    e = np.sort(e[:, e[0] != e[1]], axis=0)  # point every edge to the higher index: acyclic
    s = rng.random(e.shape[1], dtype=np.float32)
    out = get_tracks(200, e, s, backend, cc_cut=0.3, junction_cut=0.2, min_chain_length=1)
    assert all(len(t) >= 2 for ts in out.values() for t in ts)
