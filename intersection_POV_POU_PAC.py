"""Utilidad: especies con presencia en PAC, POV y POU (misma regla que sampling_effort)."""
import pandas as pd

from sampling_effort import resolve_eval_species


def intersection_pac_pov_pou() -> pd.Index:
    species = resolve_eval_species("intersection_pac_pov_pou")
    assert species is not None
    return species


if __name__ == "__main__":
    sp = intersection_pac_pov_pou()
    print(f"Intersection PAC ∩ POV ∩ POU (presence > 0): {len(sp)} species")
