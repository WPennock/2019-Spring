#Design of Flocculation Experiments

## Imports
```python
from aguaclara.play import *
import aguaclara.research.peristaltic_pump as pump
import aguaclara.research.floc_model as floc
import aguaclara.research.stock_qc as stock
import aguaclara.research.environmental_processes_analysis as epa
import aguaclara.research.procoda_parser as pro
import scipy.stats as stats
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
def v_capture(Q,D,L,angle):
    """
    """
    return ((4*Q)/(np.pi*D**2*(L/D*np.cos(angle)+np.sin(angle)))).to(u.mm/u.s)
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
v_c.to(u.mm/u.s)/Q_S.to(u.mL/u.s)
Q_S.to(u.mL/u.s)/v_c.to(u.mm/u.s)

# Residence Times
A_S = np.pi*D_S**2/4
V_S = A_S*L_S
T_S = V_S/Q_S
T_S.to(u.min)

(T_S.to(u.s)*v_c.to(u.mm/u.s))**(-1)/3
# Residence time of top of settler to turbidity meter.
T_Tube = (52*u.inch*np.pi/4*(0.25*u.inch)**2/Q_S).to(u.s)
(T_Tube/T_S).to(u.dimensionless)
# Residence time in tube after settler is about 9% of total residence time, so it will not be counted.

# Calibrate SWaT Pump
T_100mL_4_25_19 = 30*u.s
T_100mL_4_26_19 = 31*u.s

mL_rev_S = (100*u.mL/T_100mL_4_26_19)/(50.1*u.rpm)
mL_rev_S.to(u.ml/u.rev)

# SWaT Pump RPM
rpm_S = rpm_pump(Q_S,mL_rev_S)
rpm_S

v_c_Max = v_capture(max_rpm*mL_rev_S,D_S,L_S,a_S)
v_c_Max
v_c_Min = v_capture(min_rpm*mL_rev_S,D_S,L_S,a_S)
v_c_Min
v_c_now = v_capture(18*u.rpm*mL_rev_S,D_S,L_S,a_S)
v_c_now
18*u.rpm*mL_rev_S
50*u.rpm*mL_rev_S

(floc.g_straight(Q_S,0.25*u.inch)).to(1/u.s)
## SWaT Pump Calibration
```

## Experiment Design
Assuming that experiments will be conducted at a single coagulant dose for all capture velocities (0.1, 0.2, 0.3, 0.4, 0.5) mm/s, the total residence time should be at least two flocculator residence times plus triple each SWaT residence time plus 500 s to read the data.
```python
T_tot = (2*T + 3*np.sum(T_S) + 500*u.s*len(T_S)).to(u.min)
T_tot
```

## Clay Pump
```python
# Constants
C_C = 90*u.NTU # Desired clay concentration
ID_C = 2.79*u.mm # Nominal inner diameter of clay pump tubing
V_CS = 2*u.L # Volume of clay stock

# Calculations
mL_rev_nom_C = pump.vol_per_rev_3_stop(inner_diameter=ID_C)

C_CS_range = C_range(Q,C_C,mL_rev_nom_C)
C_CS_range

C_CS = 100*u.g/u.L
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
                          / floc.ratio_clay_sphere(RatioHeightDiameter))
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
Al2O3 = 0.106 # Percentage of Al2O3
SG_PSS = 1.27 # Specific gravity of stock
R_Al_Al2O3 = 0.52925 # Ratio of MW of Al to Al2O3
C_PSS = R_Al_Al2O3*Al2O3*SG_PSS*pc.density_water(293.15*u.degK) # Super stock concentration, assuming 20°C
C_PSS.to(u.g/u.L)
C_PSS = 70.6*u.g/u.L # Super stock concentration
ID_P = 1.52*u.mm
V_PS = 1*u.L

# Calculations
## Nominal capacity
mL_rev_nom_P = pump.vol_per_rev_3_stop(inner_diameter=ID_P)
mL_rev_nom_P
## Measured capacity calibration
Path_PACl = r"C:\Users\whp28\Google Drive\AGUACLARA DRIVE\AguaClara Grads\William Pennock\2019 Spring\Experiments\Data\4-25-2019\PACl Calibration\PACl_Calibration_2.xls"
PACl_Time = pro.column_of_time(Path_PACl,1)
PACl_Balance = pro.column_of_data(Path_PACl,1,7)*u.g
linreg = stats.linregress(PACl_Time,PACl_Balance)
slope, int, r_value = linreg[0:3]
Q_PACl_Cal = ((-slope*u.g/u.day)/pc.density_water(283.15*u.degK)).to(u.mL/u.min)
mL_rev_P = (Q_PACl_Cal/(50*u.rpm)).to(u.mL/u.rev)
mL_rev_P
C_CP_range = C_range(Q,C_P,mL_rev_nom_P)
C_CP_range.magnitude

C_PS = 2*u.g/u.L
V_PSS = C_PS*V_PS/C_PSS
V_PSS.to(u.mL)
Q_PS = Q_Stock(C_P,Q,C_PS)
Q_PS
rpm_PS = rpm_pump(Q_PS,mL_rev_P)
rpm_PS

T_PS = T_Stock(C_P,Q,C_PS,V_PS)
T_PS
```

## Base Pump
```python
# Compensating for PACl
n_PACl = 1.77*u.eq/u.L # Normality of PACl, eqv/L
n_NaOH = 1*u.eq/u.L # Normality of NaOH, eqv/L
eqv_PACl_m = n_PACl/C_PSS # equivalence of PACl by mass, eqv/g
eqv_PACl_m.to(u.eq/u.g)
eqv_PACl_v = C_PS*eqv_PACl_m # equpivalence of PACl by volume, eqv/L
eqv_PACl_v.to(u.eq/u.L)
V_BS = 1*u.L # Volume of base stock
MW_NaOH = 40*u.g/u.mol # Molecular weight of NaOH
m_B = MW_NaOH*V_BS*eqv_PACl_v
m_B.to(u.g)

# Acid Neutralizing Capacity
T_Alk_CaCO3 = 117*u.mg/u.L
MW_CaCO3 = 100.0869*u.g/u.mol
T_Alk = T_Alk_CaCO3/MW_CaCO3
T_Alk.to(u.eq/u.L)
# T_Alk = 2.8E-3*u.eq/u.L # eqv/L
pH = 7.508


pH_Target = 7.5 # Target pH
pH_Current = 7.9

# def M_OH(pH):
#     """This function calculates the molarity of OH from the pH.
#     Parameters
#     ----------
#     pH : float
#         pH to be inverted
#     Returns
#     -------
#     The molarity of OH (in moles per liter) of the given pH
#     Examples
#     --------
#     >>> M_OH(8.25)
#     1.778279410038923e-06 mole/liter
#     >>> M_OH(10)
#     1e-4 mole/liter
#     """
#     return 10**(pH-14)*u.mol/u.L

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
Base_Water.to(u.eq/u.L)
rpm_target = 15*u.rpm
mL_rev_nom_B = 0.21*u.mL/u.rev
# Calibrate Base pump
Path_Base = r"C:\Users\whp28\Google Drive\AGUACLARA DRIVE\AguaClara Grads\William Pennock\2019 Spring\Experiments\Data\4-25-2019\Base Calibration\Base_Calibration_2.xls"
PACl_Time = pro.column_of_time(Path_Base,1)
PACl_Balance = pro.column_of_data(Path_Base,1,7)*u.g
linreg_B = stats.linregress(PACl_Time,PACl_Balance)
slope_B, int_B, r_value_B = linreg_B[0:3]
Q_Base_Cal = ((-slope_B*u.g/u.day)/pc.density_water(283.15*u.degK)).to(u.mL/u.min)
mL_rev_B = (Q_Base_Cal/(50*u.rpm)).to(u.mL/u.rev)
mL_rev_B
N_BS2 = Base_Water*Q/(rpm_target*mL_rev_B)
N_BS2.to(u.eq/u.L)
m_B2 = MW_NaOH*V_BS*N_BS2
m_B2.to(u.g)
```
## Acid dose when pH > 7.5
```python
pH_Target = 7.5 # Target pH
pH_Current = 7.85

Acid_Water =  ANC_Current - ANC_Target
Acid_Water.to(u.meq/u.L)
N_AS = Base_Water*Q/(rpm_target*mL_rev_B)
N_AS.to(u.eq/u.L)
V_AS = 1*u.L

N_ASS = 1*u.eq/u.L
V_ASS = N_AS*V_AS/N_ASS
V_ASS.to(u.mL)
```
## Using base to solubilize PACl
It turns out that copper pipe is very insoluble at pH = 10 ([link](https://www.researchgate.net/publication/254148306_Effects_of_Changing_disinfectants_on_lead_and_copper_release/figures?lo=1)), but PACl is pretty insoluble at pH = 10 (van Benschoten and Edzwald, 1990). By adjusting the pH to 10, I aim to eliminate the film on the tubing.
```python
pH_Soluble = 10.0
ANC_Soluble = epa.ANC_closed(pH_Soluble,C_T)
ANC_Soluble.to(u.eq/u.L)
Base_Caustic = ANC_Soluble-ANC_Current
Base_Caustic.to(u.eq/u.L)
Q_B_Caustic = Base_Caustic*Q/N_BS2
Q_B_Caustic.to(u.mL/u.min)
rpm_Caustic = Q_B_Caustic/mL_rev_B
rpm_Caustic
T_Caustic = T_Stock(Base_Caustic,Q,N_BS2,V_BS)
T_Caustic

-np.log10(0.25/90)
```

##Doctest
```python
doctest.testmod(verbose=True)
```
