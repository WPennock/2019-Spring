#Design of Flocculation Experiments

## Imports
```python
from aguaclara.play import *
u.define('rev = 1 * revolutions')
import aguaclara.research.tube_sizing as ts
import aguaclara.research.floc_model as floc
import doctest
import pdb
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
## New Functions
```python
def Q_SWAT(v_c,D,L,angle):
    """This function takes the desired capture velocity, diameter of the tube settler, length of the tube settler, and its angle with respect to the horizon to give a needed flow rate through the tube settler.
    >>> from aguaclara.play import*
    >>> Q_SWAT(0.12*u.mm/u.s,1.049*u.inch,86*u.cm,60*u.deg)
    <Quantity(68.26554885933704, 'milliliter / minute')>

    """
    return ((v_c*D**2*(L/D*np.cos(angle)+np.sin(angle)))/(4/np.pi)).to(u.mL/u.min)

def Q_Stock(C, Q_plant, C_stock):
    """This function calculates the flow of the stock based on the required concentration, the flow of the plant, and the concentration of the stock.

    Parameters
    ------------
    Q_plant : float
        flow rate of the plant
    C : float
        desired concentration within the plant
    C_stock: float
        concentration of the stock

    Returns
    ----------    
    float
        Required flow rate from the stock to achieve the specified concentration at the specified flow rate, given the concentration of the stock.

    Examples
    ----------    
    >>> from aguaclara.play import*
    >>> Q_Stock(5*u.mg/u.L,1*u.L/u.s,500*u.mg/u.L)
    <Quantity(600.0, 'milliliter / minute')>
    """
    return (Q_plant*C/C_stock).to(u.mL/u.min)

def T_Stock(C, Q_plant, C_stock, V_stock):
    """ This function calculates the time it will take for a stock to run out based on the desired concentration, the flow rate of the plant, the concentration of the stock, and the volume of the stock.

    Parameters
    ------------
    Q_plant : float
        flow rate of the plant
    C : float
        desired concentration within the plant
    C_stock: float
        concentration of the stock
    V_stock: float
        volume of the stock tank

    Returns
    ----------    
    float
        Time until stock runs out

    Examples
    ----------    
    >>> from aguaclara.play import*
    >>> T_Stock(5*u.mg/u.L,1*u.L/u.s,500*u.mg/u.L,20*u.L)
    <Quantity(33.333333333333336, 'minute')>
    """
    return (V_stock/Q_Stock(C, Q_plant, C_stock)).to(u.min)

def rpm_pump(Q_stock, mL_rev):
    """This function gives the flow rate necessary for a 6 roller pump head to deliver the given flow rate for the tubing capacity given in mL/s

    Examples
    ----------    
    >>> from aguaclara.play import*
    >>> rpm_pump(10*u.mL/u.min,1*u.mL/u.rev)
    <Quantity(10.0, 'rev / minute')>
    """
    return (Q_stock/mL_rev).to(u.rev/u.min)

def C_range(Q_plant,C,mL_rev):
    """Function to determine the minimum concentration and maximum concentration (in g/L) a stock can have given the desired concentration and flow rate in the system and the capacity of the tubing (mL/rev) pumping from the stock

    Examples
    ----------    
    >>> from aguaclara.play import*
    >>> C_range(0.1*u.L/u.s,50*u.mg/u.L,0.10*u.mL/u.rev)
    <Quantity([  31.57894737 1000.        ], 'gram / liter')>
    """
    if type(C.magnitude) == np.ndarray:
      # pdb.set_trace()
      spread = np.zeros((len(C),2))
      for i in range(0,len(C)):
        spread[i,0] = (((Q_plant*C[i])/(ts.max_rpm*mL_rev)).to(u.g/u.L)).magnitude
        spread[i,1] = (((Q_plant*C[i])/(ts.min_rpm*mL_rev)).to(u.g/u.L)).magnitude
    else:
      mini = (((Q_plant*C)/(ts.max_rpm*mL_rev)).to(u.g/u.L)).magnitude
      maxi = (((Q_plant*C)/(ts.min_rpm*mL_rev)).to(u.g/u.L)).magnitude
      spread = [mini, maxi]
    return np.array(spread)*u.g/u.L
C_test = np.array([2,5,10])*u.mg/u.L    
C_range(0.1*u.L/u.s,50*u.mg/u.L,0.10*u.mL/u.rev)    
C_range(0.1*u.L/u.s,C_test,0.10*u.mL/u.rev)
```

## SWAT Pump
```python
# Constants
L_S = 86*u.cm
D_S = 1.049*u.inch
a_S = 60*u.deg

v_c = np.array([0.1,0.2,0.3,0.4,0.5])*u.mm/u.s

Q_S = Q_SWAT(v_c,D_S,L_S,a_S)
Q_S

# Residence Times
A_S = np.pi*D_S**2/4
V_S = A_S*L_S
T_S = V_S/Q_S
T_S.to(u.min)

# Need to add residence time of top of settler to turbidity meter.
```

## Experiment Design
Assuming that experiments will be conducted at a single coagulant dose for all capture velocities (0.1, 0.2, 0.3, 0.4, 0.5) mm/s, the total residence time should be at least two flocculator residence times plus double each SWaT residence time.
```python
T_tot = (2*T + 4*np.sum(T_S)).to(u.min)
T_tot
```

## Clay Pump
```python
# Constants
C_C = 50*u.NTU # Desired clay concentration
ID_C = 2.79*u.mm # Nominal inner diameter of clay pump tubing
V_CS = 4*u.L # Volume of clay stock

# Calculations
mL_rev_nom_C = ts.Q6_roller(ID_C)

C_CS_range = C_range(Q,C_C,mL_rev_nom_C)
C_CS_range

C_CS = 50*u.g/u.L
Q_CS = Q_Stock(C_C,Q,C_CS)
Q_CS
rpm_CS = rpm_pump(Q_CS,mL_rev_nom_C)
rpm_CS

T_CS = T_Stock(C_C,Q,C_CS,V_CS)
T_CS.to(u.hr)
```

## PACl Pump
```python
# Constants
Gammas = np.array([0.0,0.1,0.2,0.3,0.4,0.5])
C_P = np.zeros(len(Gammas))
for i in range(0,len(Gammas)):
  C_P[i]=0
  while np.abs(floc.gamma_coag(C_C,(C_P[i]*u.mg/u.L),floc.PACl,floc.Clay,D,floc.RATIO_HEIGHT_DIAM)-Gammas[i])>0.001:
    C_P[i]=C_P[i]+0.01
C_P = C_P*u.mg/u.L    
floc.gamma_coag(C_C,C_P,floc.PACl,floc.Clay,D,floc.RATIO_HEIGHT_DIAM)  
C_PSS = 70.6*u.g/u.L # Super stock concentration
ID_P = 1.52*u.mm
V_PS = 1*u.L

# Calculations
mL_rev_nom_P = ts.Q6_roller(ID_P)

C_CP_range = C_range(Q,C_P,mL_rev_nom_P)
C_CP_range.magnitude

C_PS = 1*u.g/u.L
Q_PS = Q_Stock(C_P,Q,C_PS)
Q_PS
rpm_PS = rpm_pump(Q_PS,mL_rev_nom_P)
rpm_PS

T_PS = T_Stock(C_P,Q,C_PS,V_PS)
T_PS
```
## Base Pump
```python
# Compensating for PACl
n_PACl = 1.77/u.L # Normality of PACl
n_NaOH = 1/u.L
eqv_PACl_m = n_PACl/C_PSS
eqv_PACl_v = C_PSS*eqv_PACl_m
eqv_PACl_v
V_BS = 1*u.L
MW_NaOH = 40*u.g
m_B = MW_NaOH*V_BS*eqv_PACl_v
m_B
```

##Doctest
```python
doctest.testmod(verbose=True)
```
