"""`--preset` resolution for dd_afpocket-sample/dd_afpocket-run: pure argparse logic, no
OpenMM/fpocket dependency, so these run everywhere the rest of the test
suite does (unlike an actual sampling run)."""
from dd_afpocket.cli import PRESETS, _resolve_preset_args, build_run_parser, build_sample_parser


def test_default_preset_matches_sample_pocket_defaults():
    args = build_sample_parser().parse_args(["p.pdb", "pocket_dir", "-o", "out"])
    assert _resolve_preset_args(args) == PRESETS["default"]


def test_quick_preset_overrides_all_bundle_keys():
    args = build_sample_parser().parse_args(["p.pdb", "pocket_dir", "-o", "out", "--preset", "quick"])
    resolved = _resolve_preset_args(args)
    assert resolved == PRESETS["quick"]
    # quick trades sampling length/replica count for speed but leaves the
    # timestep at its empirically-checked-stable value (see restraints.py).
    assert resolved["timestep_fs"] == PRESETS["default"]["timestep_fs"]
    assert resolved["sample_ns"] < PRESETS["default"]["sample_ns"]
    assert resolved["n_replicas"] < PRESETS["default"]["n_replicas"]
    # quick swaps to a cheaper GB model (measured faster than gbn2)
    assert resolved["implicit_solvent"] == "obc2"
    assert PRESETS["default"]["implicit_solvent"] == "gbn2"


def test_explicit_flag_overrides_preset_value():
    args = build_sample_parser().parse_args(
        ["p.pdb", "pocket_dir", "-o", "out", "--preset", "quick", "--n-replicas", "8"],
    )
    resolved = _resolve_preset_args(args)
    assert resolved["n_replicas"] == 8
    # every other bundled key still comes from the preset
    assert resolved["sample_ns"] == PRESETS["quick"]["sample_ns"]
    assert resolved["equil_ps"] == PRESETS["quick"]["equil_ps"]
    assert resolved["implicit_solvent"] == PRESETS["quick"]["implicit_solvent"]


def test_explicit_implicit_solvent_overrides_preset_value():
    args = build_sample_parser().parse_args(
        ["p.pdb", "pocket_dir", "-o", "out", "--preset", "quick", "--implicit-solvent", "hct"],
    )
    resolved = _resolve_preset_args(args)
    assert resolved["implicit_solvent"] == "hct"
    assert resolved["n_replicas"] == PRESETS["quick"]["n_replicas"]


def test_run_parser_has_the_same_preset_bundle():
    args = build_run_parser().parse_args(["O60674", "-o", "out", "--preset", "quick"])
    assert _resolve_preset_args(args) == PRESETS["quick"]


def test_unknown_preset_rejected():
    import pytest

    with pytest.raises(SystemExit):
        build_sample_parser().parse_args(["p.pdb", "pocket_dir", "-o", "out", "--preset", "bogus"])
