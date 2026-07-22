"""Self-contained (no server) py3Dmol overlay of every cluster's
representative structure, so the pocket-shape diversity `dd_afpocket-cluster`
produces can actually be looked at rather than only read off
`cluster_report.csv`'s RMSD/volume numbers.

Every cluster PDB already shares dd_afpocket's restrained-MD reference frame
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

from .pocket import PocketSelection, Residue, residue_labels, vert_pqr_path

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
        # per-chain "or" grouping in one selection) -- fine for dd_afpocket's
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
        title=f"dd_afpocket cluster overlay ({len(cluster_pdbs)} structure(s))", legend=legend, body=view._make_html(),
    )
    out_path.write_text(html)
    return out_path


def build_pocket_candidates_view(
    receptor_pdb: Path, selections: Sequence[PocketSelection], *,
    width: int = 900, height: int = 650, sphere_opacity: float = 0.6, label_residues: bool = True,
):
    """One py3Dmol view comparing the top-ranked candidate pockets in
    `selections` (see `pocket.top_pocket_candidates`) on a single receptor
    structure: translucent white cartoon for the whole receptor, opaque
    sticks over each pocket's lining residues, that pocket's fpocket
    alpha-sphere cloud (`pocketN_vert.pqr`, via `pocket.vert_pqr_path`) as a
    solid-colored sphere cluster depicting the cavity's actual empty-space
    volume (not just the residues around it), and one text label per lining
    residue (via `pocket.residue_labels`) -- all colored per rank from
    `_PALETTE` so the (often large) druggability gap between rank 1 and the
    rest reads visually, not just off `pocket_report.json`'s numbers.
    `sphere_opacity` defaults higher (0.6) than `build_cluster_overlay_view`'s
    cartoon opacity (0.5) since these spheres are the whole point of this
    view, not a secondary reference layer.
    """
    import py3Dmol

    receptor_pdb = Path(receptor_pdb)
    receptor_text = receptor_pdb.read_text()

    view = py3Dmol.view(width=width, height=height)
    view.addModel(receptor_text, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "white", "opacity": 0.85}})

    all_resi: List[int] = []
    model_idx = 1
    for i, sel in enumerate(selections):
        color = _PALETTE[i % len(_PALETTE)]
        by_chain = _residues_by_chain(sel.residues)
        all_resi += [r.resnum for r in sel.residues]

        for chain, resnums in by_chain.items():
            view.addStyle({"model": 0, "chain": chain, "resi": resnums}, {"stick": {"color": color, "radius": 0.28}})

        vert_pqr = vert_pqr_path(Path(sel.fpocket_out_dir), sel.fpocket_id) if sel.fpocket_out_dir else None
        if vert_pqr is not None and vert_pqr.exists():
            view.addModel(vert_pqr.read_text(), "pqr")
            view.setStyle({"model": model_idx}, {"sphere": {"color": color, "opacity": sphere_opacity}})
            model_idx += 1

        if label_residues:
            for label, (x, y, z) in residue_labels(receptor_pdb, sel.residues).values():
                view.addLabel(label, {
                    "position": {"x": x, "y": y, "z": z}, "backgroundColor": color, "backgroundOpacity": 0.8,
                    "fontColor": "white", "fontSize": 11, "showBackground": True, "borderThickness": 0.4,
                })

    if all_resi:
        view.zoomTo({"resi": all_resi})
    else:
        view.zoomTo()
    return view


def write_pocket_candidates_html(
    receptor_pdb: Path, selections: Sequence[PocketSelection], out_path: Path, *,
    width: int = 900, height: int = 650, sphere_opacity: float = 0.6, label_residues: bool = True,
) -> Path:
    """Render `build_pocket_candidates_view` to a standalone HTML file at
    `out_path` -- open it directly in a browser, no server needed."""
    out_path = Path(out_path)
    view = build_pocket_candidates_view(
        receptor_pdb, selections, width=width, height=height, sphere_opacity=sphere_opacity,
        label_residues=label_residues,
    )
    legend = " &nbsp; ".join(
        f'<span style="color:{_PALETTE[i % len(_PALETTE)]}">&#9632;</span> '
        f'Rank {sel.rank}: druggability={sel.druggability_score:.3f} ({len(sel.residues)} residue(s))'
        for i, sel in enumerate(selections)
    )
    html = _HTML_TEMPLATE.format(
        title=f"dd_afpocket pocket candidates ({len(selections)} of top-ranked)", legend=legend, body=view._make_html(),
    )
    out_path.write_text(html)
    return out_path
