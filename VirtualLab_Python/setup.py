from setuptools import setup, find_packages

setup(
    name="VirtualLab",            
    version="0.1.0",
    py_modules=["VirtualLab_Init"]
    packages=find_packages(),     
    install_requires=[            
        "numpy",
        "scipy",
        "matplotlib",
        "streamlit"
    ],
)
