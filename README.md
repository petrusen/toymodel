# Toy Model for Oxygen Isotope Fractionation

A simplified kinetic model (i.e. toy model) for phosphoryl-transfer reactions of the
different oxygen isotopologues and isotopomers of the four species. 

![Example Image](example.png)
*(output example)*

## 1. Description

The toy model code was developed to evaluate the oxygen isotope fractionation in all experimentally 
accessible chemical species, namely 2-naphthyl phosphate, water, 2-naphthol and phosphate 
(see [Bernet et al. 2026](https://chemrxiv.org/))


```
RO-PO3 + H2O     --(k1, k2)-->  RO-PO3-OH2
RO-PO3-OH2       --(k3, k4)-->  ROH2 + PO4
```

Each compound is represented by an **immutable tuple** of atomic weights. Hydrogen atoms are omitted.

| Compound     | Convention                 | Description         |
|--------------|-----------------------------|--------------------|
| RO-PO3       | `(12, 16, 16, 16, 16)`      | Organophosphate    |
| H2O          | `(16)`                      | Water              |
| RO-PO3-OH2   | `(12, 16, 16, 16, 16, 16)`  | Transition state   |
| ROH2         | `(12, 16)`                  | Alcohol            |
| PO4          | `(16, 16, 16, 16)`          | Phosphate          |

Isotopes are introduced by changing the elements of the tuple, e.g. `(12, 18, 16, 16, 16)` represents an organophosphate with an <sup>18</sup>O isotope substitution.

## 2. Dependencies

The code was tested on an Ubuntu 22.04 machine.

| Package       | Version |
|---------------|---------|
| scine-kinetx  | 3.1.0   |
| pyyaml        | 6.0.1   |
| networkx      | 2.5.1   |
| matplotlib    | 3.10.1  |
| numpy         | 2.2.4   |

### Installing `scine-kinetx`

`scine-kinetx` is not distributed via PyPI/conda and must be built from source:

```bash
# Create a dedicated environment
python3 -m venv env4toymodel
cd env4toymodel
source bin/activate

# Clone and install the C++ repository
git clone https://github.com/qcscine/kinetx.git
cd kinetx
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local -DSCINE_BUILD_PYTHON_BINDINGS=ON ..
make
make test
make install
```

This installs the Python bindings into `$HOME/.local`. Make sure this location is on your `PYTHONPATH` (and `$HOME/.local/bin` on your `PATH`, if needed):

```bash
export PYTHONPATH=$HOME/.local/lib/python3.12/site-packages:$PYTHONPATH
```

> Adjust the `python3.12` path segment to match the Python version used during the build.

After installing the Python binding of the C++ KiNetX library, it is possible to pip install the rest:

```bash
cd <path-to-the-env4toymodel>
git clone https://github.com/petrusen/toymodel.git
cd toymodel
python3 -m pip install -e .
```

## 3. Example

```bash
# Introduce input parameters
vi config.yaml  # or the text editor of your choice

# Run kinetic simulation
python3 -m toymodel
```

## Support

Should you find any problem or bug, please write a message to
[enric.petrus@eawag.ch](mailto:enric.petrus@eawag.ch) or open a new issue.
