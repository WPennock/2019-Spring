# Design Paper Calculations

## Imports
```python
from aguaclara.play import *
import aguaclara.research.floc_model as floc
import doctest
import pdb
```

## Slot Width Reynolds Number
```python
# Spacing based on H/S of  and H of 2 m
S = 40*u.cm
# Minimum velocity based on scour requirement
v = 15*u.cm/u.s
# Kinematic viscosity of water
nu = 1*u.mm**2/u.s

Re_Slot = S*v/nu
Re_Slot.to(u.dimensionless)
```
