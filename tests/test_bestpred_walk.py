from __future__ import annotations

import numpy as np
import pytest

from graph_segment import cc_and_bestpred_walk


def edges(*triples):
    e = np.array([t[:2] for t in triples]).T
    return e, np.array([t[2] for t in triples], dtype=np.float32)


def test_simple_path_is_a_cc_track():
    e, s = edges((0, 1, 0.9), (1, 2, 0.9), (2, 3, 0.9))
    assert cc_and_bestpred_walk(4, e, s) == ([[0, 1, 2, 3]], [])


def test_best_predecessor_rule():
    # chain 0-1-2-3; hit 4 has a better edge into 2 than 1 has, so the walk from 0
    # stops at 1 and 2-3 go to the walk from 4
    e, s = edges((0, 1, 0.9), (1, 2, 0.9), (2, 3, 0.9), (4, 2, 0.95))
    cc, walk = cc_and_bestpred_walk(5, e, s)
    assert cc == []
    assert sorted(walk) == [[0, 1], [4, 2, 3]]


def test_branches_keep_longest():
    # 0 -> 1 -> 2 -> 3 and 0 -> 4, both above th_add: the longer branch wins
    e, s = edges((0, 1, 0.9), (1, 2, 0.9), (2, 3, 0.9), (0, 4, 0.95))
    cc, walk = cc_and_bestpred_walk(5, e, s)
    assert walk[0] == [0, 1, 2, 3]


def test_post_filters():
    e, s = edges((0, 1, 0.9), (1, 2, 0.9), (2, 3, 0.9), (4, 5, 0.2), (5, 6, 0.2))
    assert sorted(cc_and_bestpred_walk(7, e, s)[0]) == [[0, 1, 2, 3], [4, 5, 6]]
    assert cc_and_bestpred_walk(7, e, s, min_track_score=0.5)[0] == [[0, 1, 2, 3]]
    assert cc_and_bestpred_walk(7, e, s, min_track_hits=4)[0] == [[0, 1, 2, 3]]


def test_bad_edge_index():
    with pytest.raises(ValueError, match="out of range"):
        cc_and_bestpred_walk(2, np.array([[0], [5]]), np.array([0.9], dtype=np.float32))


def test_matches_acorn_on_random_graphs():
    torch = pytest.importorskip("torch")
    pytest.importorskip("acorn")
    from torch_geometric.data import Data

    from acorn.stages.track_building.models import cc_and_walk_utils as u

    gen = torch.Generator().manual_seed(0)
    for _ in range(30):
        n, m = 300, 900
        ei = torch.randint(0, n, (2, m), generator=gen)
        ei = ei[:, ei[0] != ei[1]]
        g = Data(hit_id=torch.arange(n), hit_r=torch.rand(n, generator=gen), hit_z=torch.rand(n, generator=gen),
                 edge_index=ei, edge_scores=torch.rand(ei.shape[1], generator=gen), num_nodes=n)
        g = u.remove_cycles(g)
        for args in [(0.3, 0.4, 0.8, 0.0, 0), (0.01, 0.3, 0.8, 0.5, 4)]:
            ref = u.cc_and_bestpred_walk(g, "edge_scores", *args)
            cc, walk = cc_and_bestpred_walk(n, g.edge_index.numpy(), g.edge_scores.numpy(), *args)
            assert {frozenset(t) for t in cc} == {frozenset(t) for t in ref["cc"]}
            assert {frozenset(t) for t in walk} == {frozenset(t) for t in ref["walk"]}
