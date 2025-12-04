# -*- coding: utf-8 -*-
"""
VirtualLab_init.py

Authors: TokaLab team, 
https://github.com/TokaLab/VirtualLab
Date: 31/10/2025

This script serves as the initialization file for the TokaLab environment.
It must be run at the beginning of any session to:
1. Set up the Python paths for all project-specific modules.
2. Check that required packages are installed with the correct versions.
"""

import os
import sys
import importlib

# Name of the machine/project in use, used for paths and diagnostics
machine = "Tokalab"

# Get the main directory where this script is located
# This is used as the base for constructing all relative module paths
path_main = os.path.dirname(os.path.abspath(__file__))

# List of paths to add to Python path to allow module imports
paths_to_add = [
     "examples",
        "SimPla_Python",
        "SimPla_Python/functions",
        "SimPla_Python/tokamaks",
        "SimPla_Python/tokamaks/equilibrium",
        "SimPla_Python/tokamaks/geometry",
        "SimPla_Python/tokamaks/kinetic",
        "SynDiag_Python",
        "SynDiag_Python/diagnostics",
        f"SynDiag_Python/diagnostics/{machine}",
        f"SynDiag_Python/diagnostics/{machine}/diagnostics_data"
    ]

# Correct path separators for Unix systems (if needed)
if os.sep == '/':
    paths_to_add = paths_to_add.replace('\\', '/')

# Add paths to sys.path only if they are not already present
for relative_path in paths_to_add:
    full_path = os.path.join(path_main, relative_path)
    full_path = os.path.normpath(full_path)
    if full_path not in sys.path:
        sys.path.append(full_path)
        print(f"new added path : {full_path}")
        
"""
Environment check — verifies that required packages are installed
and that their versions match the ones needed for TokaLab to work correctly.
"""

# Dictionary of required packages and their exact versions
required_packages = {
    "numpy": "1.26.4",
    "scipy": "1.15.3",
    "matplotlib": "3.9.2"}

# List to collect errors or warnings
errors = []

# Check each package
for pkg, required_version in required_packages.items():
    try:
        # Dynamically import the package
        module = importlib.import_module(pkg)
        
        # Get the installed version of the package
        installed_version = getattr(module, "__version__", None)
        
        if installed_version is None:
            errors.append(f"⚠️  {pkg}: unable to determine version")
        elif installed_version != required_version:
            errors.append(
                f"⚠️  {pkg}: installed version {installed_version}, "
                f"required version {required_version}"
            )
    except ImportError:
        errors.append(f"❌  {pkg}: not installed")

# Print the results
if errors:
    print("\n=== Environment check — problems found ===")
    for e in errors:
        print(e)
    print("\n➡️ Please consider that you have some packages with different versions." 
          "If you encounter some errors, try to install the correct version\n")

else:
    print("✅ All required packages are installed with correct versions.")



