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
# Nominal plane jet outlet dimension
b = 0.37*S
# Minimum velocity based on scour requirement
v = 15*u.cm/u.s
# Nominal initial velocity of plane jet
U_0 = v/0.37
# Kinematic viscosity of water
nu = 1*u.mm**2/u.s

Re_Plane_Jet = b*U_0/nu
Re_Plane_Jet.to(u.dimensionless)
```
