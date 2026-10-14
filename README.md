# CO2 Chedin Potential

Python implementation of Chedin's spectroscopic potential for the CO2 molecule in quasi-normal symmetry coordinates.

This repository provides a compact, high-performance implementation of a CO2 intramolecular potential energy model suitable for electronic-structure-inspired force evaluations and minimal molecular simulations. The project combines a Python/ASE-based implementation with a compiled C extension for faster energy and force calculations.

## Overview

The model follows the CO2 potential developed by Chedin and related work on quasi-normal symmetry coordinates. The implementation evaluates the potential energy and its forces as a function of the C-O bond distances and bend angle, using the polynomial expansion in symmetry coordinates used for spectroscopic modeling.

The code is structured as follows:

- `chedin_pot/` — Python package and compiled C extension for performance-critical energy/force evaluation.
- `c_code/` — reference or legacy C implementation files.
- `References/` — papers and background literature used for the model.
- `co2.ipynb` — notebook demonstrating usage, benchmarking, and a short molecular dynamics example.

## Features

- CO2 potential energy model based on Chedin's quasi-normal symmetry coordinates.
- Force evaluation for ASE-compatible calculations.
- High-performance C extension built with NumPy and setuptools.
- Example notebook for validation and benchmarking.
- Support for periodic boundary conditions in the compiled extension implementation.

## Project layout

```text
.
├── README.md
├── co2.ipynb
├── c_code/
│   ├── co2.ccl4.in
│   ├── co2.in.nm
│   └── intrapot.c
├── chedin_pot/
│   ├── __init__.py
│   ├── chedin.c
│   ├── pyproject.toml
│   └── setup.py
├── References/
│   ├── Chedin-1979-JMolSpectr-76-430.pdf
│   ├── Chedin-1984-JMolSpectr-107-333.pdf
│   ├── 2003-NPCS-6-502.pdf
│   ├── 2005-PCCP-7-1205.pdf
│   └── phys01.pdf
└── ...
```

## Installation

This project requires Python and NumPy. The compiled extension is built using setuptools.

From the repository root:

```bash
cd chedin_pot
python -m pip install -e .
```

This will build the `ChedinPot_calculator` extension in editable mode.

## Usage

After installation, the compiled C extension can be imported directly or through the package wrapper:

```python
from chedin_pot import compute
```

A minimal ASE-based example is shown below:

```python
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from chedin_pot import compute

class Chedin(Calculator):
    implemented_properties = ['energy', 'forces']

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        positions = np.ascontiguousarray(self.atoms.get_positions(), dtype=np.float64)
        cell = np.ascontiguousarray(self.atoms.get_cell(), dtype=np.float64)
        forces = np.ascontiguousarray(np.zeros((len(self.atoms), 3)), dtype=np.float64)

        energy = compute(positions, cell, forces)
        self.results['energy'] = energy * 0.010364
        self.results['forces'] = forces * 0.010364

atoms = Atoms(
    'CO2',
    positions=[[0.0, 0.0, 0.0], [1.161, 0.0, 0.0], [-1.161, 0.0, 0.0]],
    cell=np.eye(3) * 10.0,
    pbc=True,
)

atoms.calc = Chedin()
print(atoms.get_potential_energy())
print(atoms.get_forces())
```

A complete working notebook is available in `co2.ipynb`, including:

- a Python implementation of the same potential,
- a compiled C version,
- a force comparison between the two,
- and a short MD example with ASE.

## Notes on the model

The implementation uses quasi-normal symmetry coordinates for a linear triatomic molecule, with a polynomial expansion in terms of the symmetric and antisymmetric stretching/bending coordinates. The energy expression includes quadratic, cubic, quartic, quintic, and sextic contributions, following the form used in Chedin's spectroscopic potential model.

This makes the model suitable for exploring CO2 structural energetics and force evaluations in a compact, analytic form.

## References

The repository includes reference articles in the `References/` directory, including works by:

- Chedin, A.; Teffo, J.-L. J. Mol. Spectr. 107 (1984) 333
- Simons, G.; Parr, R.G.; Finlan, J.M. J. Chem. Phys. 59 (1973) 3229
- Carney, G.D.; Sprandel, L.L.; Kern, C.W. Adv. Chem. Phys. 37 (1977) 305

## License

This project does not currently declare an explicit license file. Please check the repository contents and any accompanying project metadata before using it in a production or publication context.

## Contributing

Contributions, fixes, and extensions are welcome. If you improve the implementation, add tests, or document new use cases, please open a pull request with a clear explanation of the change.
