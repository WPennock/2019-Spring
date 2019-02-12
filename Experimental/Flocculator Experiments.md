#Design of Flocculation Experiments

## Imports
```python
from aguaclara.play import *
u.define('rev = 1 * revolutions')
import aguaclara.research.tube_sizing as ts
import doctest
```

## System Constants
```python
Q = 100*u.mL/u.s
L = 56.35*u.m
D = 1.25*u.inch
A = np.pi*D**2/4
V = L*A
T = V/Q
T.to(u.min)
# The nominal residence time of the flocculator is 7.5 minutes.
```

## SWAT Pump
```python
# Constants
L_S = 86*u.cm
D_S = 1.049*u.inch
a_S = 60*u.deg

v_c = np.array([0.1,0.2,0.3,0.4,0.5])*u.mm/u.s

def Q_SWAT(v_c,D,L,angle):
    """This function takes the desired capture velocity, diameter of the tube settler, length of the tube settler, and its angle with respect to the horizon to give a needed flow rate through the tube settler.
    >>> from aguaclara.play import*
    >>> Q_SWAT(0.12*u.mm/u.s,1.049*u.inch,86*u.cm,60*u.deg)
    <Quantity(68.26554885933704, 'milliliter / minute')>

    """
    return ((v_c*D**2*(L/D*np.cos(angle)+np.sin(angle)))/(4/np.pi)).to(u.mL/u.min)

Q_S = Q_SWAT(0.12*u.mm/u.s,D_S,L_S,a_S)
Q_S

# Residence Times
A_S = np.pi*D_S**2/4
V_S = A_S*L_S
T_S = V_S/Q_S
T_S.to(u.min)
```

## Clay Pump
```python
C_C = 50*u.NTU
ID_C = 2.79*u.mm

C_CS_Max = C_C*Q/(ts.Q6_roller(ID_C)*ts.min_rpm)
C_CS_Max.to(u.g/u.L)

C_CS_Min = C_C*Q/(ts.Q6_roller(ID_C)*ts.max_rpm)
C_CS_Min.to(u.g/u.L)

C_CS = 100*u.g/u.L
T_CS = ts.T_stock(Q,C_C,)
```

## Example code from Fluoride Auto
```python
#Given: flow rate of system, pump speed from stock, concentration of stock
#Find: system concentration

#Q_sys * C_sys = Q_stock * C_stock

pump_speed = 30*(u.rpm)
orange_yellow = 0.019*(u.milliliter/u.revolutions)
oy_flowrate = orange_yellow.to(u.milliliter/u.revolutions)*(pump_speed).to(u.revolutions/u.s)

Q_sys = 0.7601 * (u.mL/u.s)
C_stock = 2400 * (u.mg/u.L)
C_sys = 10 * (u.mg/u.L)


Q_stock = (Q_sys * C_sys)/C_stock
Q_stock_rpm = Q_stock *
```
##Doctest
```python
doctest.testmod(verbose=True)
```
