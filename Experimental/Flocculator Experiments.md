#Design of Flocculation Experiments

## Imports
```python
from aguaclara.play import *
u.define('rev = 1 * revolutions')
import aguaclara.research.peristaltic_pump as pump
import aguaclara.research.floc_model as floc
import aguaclara.research.stock_qc as stock
import aguaclara.research.environmental_processes_analysis as epa
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
max_rpm = 100*u.rev/u.min
min_rpm = 3*u.rev/u.min

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
        spread[i,0] = (((Q_plant*C[i])/(max_rpm*mL_rev)).to(u.g/u.L)).magnitude
        spread[i,1] = (((Q_plant*C[i])/(min_rpm*mL_rev)).to(u.g/u.L)).magnitude
    else:
      mini = (((Q_plant*C)/(max_rpm*mL_rev)).to(u.g/u.L)).magnitude
      maxi = (((Q_plant*C)/(min_rpm*mL_rev)).to(u.g/u.L)).magnitude
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
mL_rev_nom_C = pump.vol_per_rev_3_stop(inner_diameter=ID_C)

C_CS_range = C_range(Q,C_C,mL_rev_nom_C)
C_CS_range

C_CS = 50*u.g/u.L
m_CS = V_CS*C_CS
Q_CS = Q_Stock(C_C,Q,C_CS)
Q_CS
rpm_CS = rpm_pump(Q_CS,mL_rev_nom_C)
rpm_CS

T_CS = T_Stock(C_C,Q,C_CS,V_CS)
T_CS.to(u.hr)
```

## PACl Pump
```python
MOLEC_WEIGHT_ALUMINUM = 0.027*u.kg/u.mol

def conc_precipitate(ConcAluminum, coag):
    """Return coagulant precipitate concentration given aluminum dose.
    Note that conc_precipitate returns a value that varies from the equivalent
    MathCAD function beginning at the third decimal place. The majority of
    functions below this point in the file ultimately call on conc_precipitate
    at some point, and will not return the same value as their equivalent
    function in MathCAD. This is known.
    """
    MOLEC_WEIGHT_ALUMINUM = 0.027*u.kg/u.mol
    return (ConcAluminum * coag.PrecipMolecWeight) / (MOLEC_WEIGHT_ALUMINUM * coag.PrecipAluminumMPM)
conc_precipitate(1*u.mg/u.L,floc.PACl)

def frac_vol_floc_initial(ConcAluminum, ConcClay, coag, material):
        """Return the fraction of flocs initially present."""
        return ((conc_precipitate(ConcAluminum, coag)/coag.PrecipDensity)
                + (ConcClay / material.Density))

frac_vol_floc_initial(1*u.mg/u.L,1*u.mg/u.L,floc.PACl,floc.Clay)

def ratio_area_clay_total(ConcClay, material, DiamTube, RatioHeightDiameter):
    """Return the surface area of clay normalized by total surface area.
    Total surface area is a combination of clay and reactor wall
    surface areas. This function is used to estimate how much coagulant
    actually goes to the clay.
    """
    return (1
            / (1
               + (2 * material.Diameter
                  / (3 * DiamTube * floc.ratio_clay_sphere(RatioHeightDiameter)
                     * (ConcClay / material.Density)
                     )
                  )
               )
)

def gamma_coag(ConcClay, ConcAluminum, coag, material,
               DiamTube, RatioHeightDiameter):
    """Return the coverage of clay with nanoglobs.
    This function accounts for loss to the tube flocculator walls
    and a poisson distribution on the clay given random hits by the
    nanoglobs. The poisson distribution results in the coverage only
    gradually approaching full coverage as coagulant dose increases.
    """
    # pdb.set_trace()
    return (1 - np.exp((
                       (-frac_vol_floc_initial(ConcAluminum, 0*u.kg/u.m**3, coag, material)
                         * material.Diameter)
                        / (frac_vol_floc_initial(0*u.kg/u.m**3, ConcClay, coag, material)
                           * coag.Diameter))
                       * (1 / np.pi)
                       * (ratio_area_clay_total(ConcClay, material,
                                                DiamTube, RatioHeightDiameter)
                          / ratio_clay_sphere(RatioHeightDiameter))
))

# Constants
Gammas = np.array([0.0,0.1,0.2,0.3,0.4,0.5])
C_P = np.zeros(len(Gammas))
for i in range(0,len(Gammas)):
  C_P[i]=0
  while np.abs(gamma_coag(C_C,(C_P[i]*u.mg/u.L),floc.PACl,floc.Clay,D,floc.RATIO_HEIGHT_DIAM)-Gammas[i])>0.001:
    C_P[i]=C_P[i]+0.01
C_P = C_P*u.mg/u.L    
gamma_coag(C_C,C_P,floc.PACl,floc.Clay,D,floc.RATIO_HEIGHT_DIAM)  
C_PSS = 70.6*u.g/u.L # Super stock concentration
ID_P = 1.52*u.mm
V_PS = 1*u.L

# Calculations
mL_rev_nom_P = pump.vol_per_rev_3_stop(inner_diameter=ID_P)

C_CP_range = C_range(Q,C_P,mL_rev_nom_P)
C_CP_range.magnitude

C_PS = 1*u.g/u.L
V_PSS = C_PS*V_PS/C_PSS
V_PSS.to(u.mL)
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
n_PACl = 1.77*u.equivalents/u.L # Normality of PACl, eqv/L
n_NaOH = 1*u.eq/u.L # Normality of NaOH, eqv/L
eqv_PACl_m = n_PACl/C_PSS # equivalence of PACl by mass, eqv/g
eqv_PACl_v = C_PSS*eqv_PACl_m # equpivalence of PACl by volume, eqv/L
eqv_PACl_v
V_BS = 1*u.L # Volume of base stock
MW_NaOH = 40*u.g/u.mol # Molecular weight of NaOH
m_B = MW_NaOH*V_BS*eqv_PACl_v
m_B.to(u.g)

# Acid Neutralizing Capacity
T_Alk = 2.8E-3*u.mol/u.L # eqv/L
pH = 7.67
pH_Target = 7.5 # Target pH
pH_Current = 7.1

def M_OH(pH):
    """This function calculates the molarity of OH from the pH.
    Parameters
    ----------
    pH : float
        pH to be inverted
    Returns
    -------
    The molarity of OH (in moles per liter) of the given pH
    Examples
    --------
    >>> M_OH(8.25)
    1.778279410038923e-06 mole/liter
    >>> M_OH(10)
    1e-4 mole/liter
    """
    return 10**(pH-14)*u.mol/u.L

def Total_Carbonates(pH, Total_Alkalinity):
    """Total carbonates (C_T) calculated from pH and total alkalinity.
    Parameters
    ----------
    pH : float
        pH of the system
    Total_Alkalinity: float
        total alkalinity of the system

    Returns
    -------
    total carbonates in the system (mole/L)

    Examples
    --------
    >>> Total_Carbonates(8.25, 1*u.mol/u.L)
    0.9968913136948984 equivalents/liter
    >>> Total_Carbonates(10, 1*u.mol/u.L)
    1.359830978063728 equivalents/liter
    """
    return (Total_Alkalinity + epa.invpH(pH) - epa.Kw /
        epa.invpH(pH)) / (2*epa.alpha2_carbonate(pH) + epa.alpha1_carbonate(pH))

C_T = Total_Carbonates(pH,T_Alk).to(u.meq/u.L)        
ANC_Current = epa.ANC_closed(pH_Current,C_T)
ANC_Current.to(u.meq/u.L)
ANC_Target = epa.ANC_closed(pH_Target,C_T)
ANC_Target.to(u.meq/u.L)

Base_Water = ANC_Target - ANC_Current
rpm_target = 15*u.rpm
mL_rev_nom_B = 0.21*u.mL/u.rev

N_BS = Base_Water*Q/(rpm_target*mL_rev_nom_B)
N_BS.to(u.eq/u.L)
```

##Doctest
```python
doctest.testmod(verbose=True)
```
