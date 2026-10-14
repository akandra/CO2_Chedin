import os
import numpy as np
from setuptools import setup, Extension

# 1. Define the C extension module
ase_c_extension = Extension(
    name='ChedinPot_calculator',               # The name you will use to import it in Python
    sources=['chedin.c'],              # Your C source code files
    include_dirs=[np.get_include()],      # Automatically inject paths for NumPy C headers
    extra_compile_args=['-O3', '-march=native'], # High-performance optimization flags for Linux
)

# 2. Execute the setup sequence
setup(
    name='ChedinPot_calculator',
    version='0.1.0',
    description='High-performance C extension for ASE Chedin"s energy and force',
    ext_modules=[ase_c_extension],
)
