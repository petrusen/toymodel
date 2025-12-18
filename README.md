# Toy model - Kinetics of oxygen isotope effects

## 1. Description

Simplistic model with many details removed so that it can be used to explain a mechanism concisely.
The simplified reaction is:

RO-PO3 + H2O --> RO-PO3-OH2 (k1, k2)
RO-PO3-OH2   --> ROH2 + PO4 (k3, k4)

Each compound is represented by a --unmutable-- tupple with the atomic weight. Hydrogen atoms
are omitted.

| Compounds       | Convention               | Description
|-----------------|--------------------------|----------------------
| RO-PO3          | (12, 16, 16, 16, 16)     | Organophosphate.
| H2O             | (16)                     | Water.
| RO-PO3-OH2      | (12, 16, 16, 16, 16, 16) | Transition state.
| ROH2            | (12, 16)                 | Alcohol.
| PO4             | (16, 16, 16, 16)         | Phosphate.

Notice that to isotops are introduced by changing the elements of the tupple, e.g. (12, 18, 16, 16, 16).

## 2. Dependencies

Below the dependencies for running the toy model are collected.

| Package         | Version |
|-----------------|---------|
| python          | 3.12.x  |
| scine-kinetx    | 3.1.0   |
| pyyaml          | 6.0.1   |
| network         | 2.5.1   |
| matplotlib      | 3.10.1  |
| numpy           | 2.2.4   |

## 3. Example

```bash

# Introduce input parameters
vi config.yaml  # or the text editor of your choice

# Run kinetic simulation
python main.py


## Support 
Should you find any problem or bug, please write a message to
[enric.petrus@eawag.ch](enric.petrus@eawag.ch) or open an new issue.
