"""Self-contained (no server) py3Dmol overlay of every cluster's
representative structure, so the pocket-shape diversity `dd_af-cluster`
produces can actually be looked at rather than only read off
`cluster_report.csv`'s RMSD/volume numbers.

Every cluster PDB already shares dd_af's restrained-MD reference frame
(`restraints.py` holds everything outside the pocket neighborhood fixed to
the pre-MD structure's own coordinates), so no additional structural
alignment step is needed before overlaying them -- unlike a typical
multi-structure comparison, where each structure would first need superposing
onto a common frame.

`py3Dmol.view._make_html()` returns a `<div>` + `<script>` fragment (not a
full document) that loads 3Dmol.js from a CDN
(https://cdn.jsdelivr.net/npm/3dmol) at view time -- opening the written
file requires internet access, the same as py3Dmol's own Jupyter usage.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence

from .pocket import Residue

# matplotlib's "tab10" palette (hex), cycled if there are more clusters than
# colors -- chosen for being visually distinct at a glance, not for any
# domain meaning.
_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

_HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>{title}</title></head>
<body>
<h3>{title}</h3>
<p>{legend}</p>
{body}
</body>
</html>
"""


def _residues_by_chain(pocket_residues: Sequence[Residue]) -> Dict[str, List[int]]:
    by_chain: Dict[str, List[int]] = {}
    for r in pocket_residues:
        by_chain.setdefault(r.chain, []).append(r.resnum)
    return by_chain


def build_cluster_overlay_view(
    cluster_pdbs: Sequence[Path], pocket_residues: Sequence[Residue] = (), *, width: int = 900, height: int = 650,
):
    """One py3Dmol view with every `cluster_pdbs` entry loaded as its own
    model: translucent cartoon for the whole (frozen) receptor, opaque
    sticks over the pocket-lining residues (`pocket_residues`) -- so the
    part that actually varies between clusters is the part that's easiest
    to see. Each cluster gets one color from `_PALETTE`, used for both its
    cartoon and its pocket sticks.
    """
    import py3Dmol

    view = py3Dmol.view(width=width, height=height)
    by_chain = _residues_by_chain(pocket_residues)

    for i, pdb_path in enumerate(cluster_pdbs):
        color = _PALETTE[i % len(_PALETTE)]
        view.addModel(Path(pdb_path).read_text(), "pdb")
        view.setStyle({"model": i}, {"cartoon": {"color": color, "opacity": 0.5}})
        for chain, resnums in by_chain.items():
            view.addStyle(
                {"model": i, "chain": chain, "resi": resnums},
                {"stick": {"color": color, "radius": 0.2}},
            )

    if by_chain:
        # Zoom to the union of pocket residues across every model, since
        # that's the region the overlay exists to compare -- the rest of
        # the receptor is identical (frozen) across clusters by
        # construction. Chain is omitted from the selector (3Dmol.js has no
        # per-chain "or" grouping in one selection) -- fine for dd_af's
        # normal single-chain AFDB inputs; a multi-chain pocket zooms to the
        # union of residue numbers across chains instead of each chain's
        # own residues, a minor imprecision only when resnums also collide
        # across chains.
        all_resnums = sorted({resnum for resnums in by_chain.values() for resnum in resnums})
        view.zoomTo({"resi": all_resnums})
    else:
        view.zoomTo()
    return view


def write_cluster_overlay_html(
    cluster_pdbs: Sequence[Path], pocket_residues: Sequence[Residue], out_path: Path, *,
    width: int = 900, height: int = 650,
) -> Path:
    """Render `build_cluster_overlay_view` to a standalone HTML file at
    `out_path` -- open it directly in a browser, no server needed."""
    out_path = Path(out_path)
    view = build_cluster_overlay_view(cluster_pdbs, pocket_residues, width=width, height=height)
    legend = ", ".join(
        f'<span style="color:{_PALETTE[i % len(_PALETTE)]}">&#9632;</span> {Path(p).stem}'
        for i, p in enumerate(cluster_pdbs)
    )
    html = _HTML_TEMPLATE.format(
        title=f"dd_af cluster overlay ({len(cluster_pdbs)} structure(s))", legend=legend, body=view._make_html(),
    )
    out_path.write_text(html)
    return out_path
