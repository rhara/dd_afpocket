"""Fixed-column PDB line builders shared across tests (same columns as the
real PDB format used elsewhere, matching dd_docking's test style)."""


def atom_line(serial, name, resname, chain, resnum, x, y, z, *, altloc=" ",
              icode=" ", element=None, rec="ATOM  ", bfactor=0.0):
    if element is None:
        element = name.strip()[0]
    name_field = f" {name:<3}" if len(name) < 4 else name
    return (f"{rec}{serial:>5} {name_field}{altloc}{resname:>3} {chain}{resnum:>4}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{bfactor:6.2f}          {element:>2}")


def ter_line(serial, resname, chain, resnum):
    return f"TER   {serial:>5}      {resname:>3} {chain}{resnum:>4}"
