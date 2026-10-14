"""Compatibility wrapper exposing the C extension as chedin_ase_calc.compute."""

from ChedinPot_calculator import compute

__all__ = ["compute"]
