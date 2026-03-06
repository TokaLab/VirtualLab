from setuptools import setup, find_packages

setup(
    name="VirtualLab",
    version="0.1.0",
    packages=find_packages(),  # trova automaticamente tutti i package con __init__.py
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "streamlit"
    ],
    py_modules=["VirtualLab_init"],  # se vuoi includere il singolo file VirtualLab_init.py
)
