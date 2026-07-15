import numpy as np
import pandas as pd
import pytest

from dd_af.cluster import cluster_frames, pick_medoids

# Two well-separated 2-point "clusters" encoded directly as a distance
# matrix (frames 0,1 near each other; frames 2,3 near each other; the two
# groups far apart) -- avoids needing a real mdtraj Trajectory for the
# pure clustering-logic tests.
_RMSD_MATRIX = np.array([
    [0.0, 0.5, 8.0, 8.2],
    [0.5, 0.0, 8.1, 8.3],
    [8.0, 8.1, 0.0, 0.3],
    [8.2, 8.3, 0.3, 0.0],
])


def test_cluster_frames_separates_well_separated_groups():
    labels = cluster_frames(_RMSD_MATRIX, n_clusters=2)
    assert labels[0] == labels[1]
    assert labels[2] == labels[3]
    assert labels[0] != labels[2]


def test_cluster_frames_caps_n_clusters_to_frame_count(capsys):
    labels = cluster_frames(_RMSD_MATRIX, n_clusters=100)
    assert len(set(labels.tolist())) == 4  # capped to n_frames
    captured = capsys.readouterr()
    assert "capping at 4" in captured.out


def test_pick_medoids_selects_frame_with_min_mean_intra_cluster_distance():
    # Cluster 0: frames 0,1 (distances 0.5 apart) -- either could tie, but
    # a 3rd nearly-identical member (frame 4) makes frame 0 the clear
    # medoid (closest to both 1 and 4).
    matrix = np.array([
        [0.0, 0.5, 9.0, 0.1],
        [0.5, 0.0, 9.1, 0.6],
        [9.0, 9.1, 0.0, 9.0],
        [0.1, 0.6, 9.0, 0.0],
    ])
    labels = np.array([0, 0, 1, 0])
    medoids = pick_medoids(matrix, labels)
    assert medoids[1] == 2  # only member of cluster 1
    assert medoids[0] == 0  # frame 0 has the lowest mean distance to {1, 3}


def test_pick_medoids_singleton_cluster_is_its_own_medoid():
    matrix = np.zeros((3, 3))
    labels = np.array([0, 1, 1])
    medoids = pick_medoids(matrix, labels)
    assert medoids[0] == 0
    assert medoids[1] in (1, 2)


def test_cluster_report_row_frame_counts_sum_to_total(tmp_path):
    # Exercises the pure bookkeeping in isolation (no real trajectory):
    # given cluster labels for N frames, the member counts summed across
    # clusters must equal N -- the same invariant
    # `cluster_pocket_trajectory` asserts on real data.
    labels = cluster_frames(_RMSD_MATRIX, n_clusters=2)
    counts = pd.Series(labels).value_counts()
    assert counts.sum() == _RMSD_MATRIX.shape[0]
