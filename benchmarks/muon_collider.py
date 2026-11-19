"""Benchmark every graph_segment backend on scored muon-collider graphs: physics and time.

Needs acorn (https://gitlab.cern.ch/gnn4itkteam/acorn) to read the graphs and match tracks.

    python benchmarks/muon_collider.py <scored dir>/valset -o out_dir -j 64
    python benchmarks/muon_collider.py --report out_dir      # table from out_dir/raw.pkl

Per event and config: get_tracks wall time (1 thread; loading, pointing the edges outward
and matching excluded), then acorn matching (ATLAS and two_way, 50% of hits, tracks of
>= 5 hits), with and without the built-in filter (min_track_score 0.5, min_track_hits 7).
"""

from __future__ import annotations

import argparse
import glob
import os
import pickle
import time
from multiprocessing import Pool

import numpy as np

CONFIGS = {  # name: (backend, cuts)
    "cc@0.01": ("cc", {"cc_cut": 0.01}),
    "cc@0.9": ("cc", {"cc_cut": 0.9}),
    "cc_walk": ("cc_walk", {"cc_cut": 0.01, "th_min": 0.3, "th_add": 0.8}),
    "cc_bestpred_walk": ("cc_bestpred_walk", {"cc_cut": 0.01, "th_min": 0.3, "th_add": 0.8}),
    "cc_junction_removal": ("cc_junction_removal",
                            {"cc_cut": 0.7, "junction_cut": 0.5, "min_chain_length": 3}),
}
FILTER = {"min_track_score": 0.5, "min_track_hits": 7}


def counts(g, trks, style):
    from acorn.stages.track_building import utils as tbu
    from acorn.stages.track_building.models.cc_and_walk_utils import add_track_labels

    g = g.clone()
    hid = g.hit_id.numpy()
    add_track_labels(g, {m: [hid[t].tolist() for t in ts] for m, ts in trks.items()})
    d = tbu.evaluate_labelled_graph(g, sel_conf={"pt": [0, float("inf")]}, matching_fraction=0.5,
                                    matching_style=style, min_track_length=5)
    real = d[d.particle_id != 0]  # hits on no true track form a fake "particle 0"
    parts = real[real.is_reconstructable]
    reco = parts[(parts.track_id != -1) & parts.is_reconstructed & parts.is_matchable]
    trk = d[(d.track_id != -1) & d.is_matchable]
    mt = real[(real.track_id != -1) & real.is_matched & real.is_matchable]
    return np.array([len(parts.drop_duplicates(["event_id", "particle_id"])),
                     len(reco.drop_duplicates(["event_id", "particle_id"])), len(reco),
                     len(trk.drop_duplicates(["event_id", "track_id"])),
                     len(mt.drop_duplicates(["event_id", "track_id"]))])


def run_event(f):
    import torch

    torch.set_num_threads(1)
    from acorn.stages.track_building.models.cc_and_walk_utils import remove_cycles
    from acorn.utils.loading_utils import add_variable_name_prefix_in_pyg, load_pyg

    from graph_segment import get_tracks

    g = add_variable_name_prefix_in_pyg(load_pyg(f.removesuffix(".gz"), weights_only=False))
    g.edge_scores = g.edge_scores.float()
    g = remove_cycles(g)
    n, ei, s = g.hit_id.shape[0], g.edge_index.numpy(), g.edge_scores.numpy()
    out = {"event": os.path.basename(f), "n_hits": n, "n_edges": ei.shape[1]}
    for name, (backend, cuts) in CONFIGS.items():
        res = {}
        for tag, extra in [("raw", {}), ("filt", FILTER)]:
            t0 = time.perf_counter()
            trks = get_tracks(n, ei, s, backend, **cuts, **extra)
            res[f"time_{tag}"] = time.perf_counter() - t0
            res[tag] = np.concatenate([counts(g, trks, "ATLAS"), counts(g, trks, "two_way")])
        out[name] = res
    return out


def report(rows, fh=None):
    p = lambda *a: print(*a, file=fh)  # noqa: E731
    p(f"{len(rows)} events; filt = {FILTER}\n")
    p(f"{'config':20s} {'filter':6s} {'eff':>6s} {'fake':>6s} {'dup':>6s} {'eff2way':>7s} "
      f"{'trk/evt':>7s} | {'t_mean_ms':>9s} {'t_med_ms':>8s} {'t_p99_ms':>8s} {'t_max_ms':>8s}")
    for n in CONFIGS:
        for tag in ["raw", "filt"]:
            c = sum(r[n][tag] for r in rows)
            t = 1e3 * np.array([r[n][f"time_{tag}"] for r in rows])
            p(f"{n:20s} {tag:6s} {c[1] / c[0]:6.3f} {1 - c[4] / c[3]:6.3f} "
              f"{(c[2] - c[1]) / max(c[1], 1):6.3f} {c[6] / c[5]:7.3f} {c[3] / len(rows):7.2f} | "
              f"{t.mean():9.2f} {np.median(t):8.2f} {np.percentile(t, 99):8.2f} {t.max():8.2f}")
    p("\ncuts: " + "; ".join(f"{n} {c}" for n, (_, c) in CONFIGS.items()))
    p("time = graph_segment.get_tracks wall time per event, 1 thread")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("input_dir", nargs="?")
    ap.add_argument("-o", "--out")
    ap.add_argument("-n", type=int, default=None)
    ap.add_argument("-j", type=int, default=32)
    ap.add_argument("--report", help="print the table from <dir>/raw.pkl")
    a = ap.parse_args()
    if a.report:
        with open(f"{a.report}/raw.pkl", "rb") as fh:
            report(pickle.load(fh))
        raise SystemExit
    os.makedirs(a.out, exist_ok=True)
    files = sorted(glob.glob(f"{a.input_dir}/*.pyg*"))[: a.n]
    with Pool(a.j, maxtasksperchild=20) as pool:
        rows = list(pool.imap_unordered(run_event, files))
    with open(f"{a.out}/raw.pkl", "wb") as fh:
        pickle.dump(rows, fh)
    with open(f"{a.out}/table.txt", "w") as fh:
        report(rows, fh)
    with open(f"{a.out}/table.txt") as fh:
        print(fh.read())
