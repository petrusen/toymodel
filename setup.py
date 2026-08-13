from setuptools import setup, find_packages
import re

with open("src/toymodel/__init__.py") as f:
    version = re.search(r'__version__ = "(.*?)"', f.read()).group(1)

setup(
    name="toymodel",
    version=version,
    description="A toy model to run kinetic simulations regarding the oxygen \
    isotope fractionation in phosphoryl transfer reactions",
    author="Enric Petrus",
    author_email=["enric.petrus@eawag.ch"],
    # Tells setuptools to look inside the 'src' directory
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    
    # Python version and dependencies
    python_requires=">=3.8",
    install_requires=[
        "scine-kinetx==3.1.0",
        "pyyaml>=6.0",
        "matplotlib>=3.5",
        "networkx>=2.5.1",
        "numpy>=2.2.4",
        "setuptools"
    ],
)
