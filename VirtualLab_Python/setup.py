from setuptools import setup, find_packages

setup(
    name="tokalab",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "streamlit",  # se lo usi
        # altre dipendenze
    ],
)