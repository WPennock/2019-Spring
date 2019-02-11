#Design of Flocculation Experiments

##Imports
```python
from aguaclara.play import *
u.define('rev = 1 * revolutions')
import aguaclara.research.tube_sizing as ts
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
