from __future__ import unicode_literals
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import matplotlib
import matplotlib.pyplot as plt
import pint 
ureg = pint.UnitRegistry(system='mks')
#ureg.load_definitions('C:\Users\William\my_def.txt') 
ureg.define('NTU = 100/68 * (mg/L)')

# Import Data
#name = raw_input('Name of File (excluding .txt): ')

ID,pH,HLf,AT,WT,ACT,ACD,Eff,Inf,pC,D,Vc,t,Qt,HLt,EDR,nu,Gam,phi,PV,PI = np.loadtxt('Agg.txt',unpack=True,skiprows=0)

# Assign Units

HLf = HLf * ureg.cm
AT = AT * ureg.degK
WT = WT * ureg.degK
ACT = ACT * ureg.degK
ACD = ACD * ureg.cm
Eff = Eff * ureg.NTU
Inf = Inf * ureg.NTU
D = D * ureg.mg / ureg.L
Vc = Vc * ureg.mm / ureg.s
t = t * ureg.s
Qt = Qt * ureg.mL / ureg.s
HLt = HLt * ureg.cm
EDR = EDR * ureg.mW / ureg.kg
nu = nu * ureg.mm**2 / ureg.s

# Calculate Composite Parameter

D_Cl = 7 * ureg.um
D_Co = 90 * ureg.nm
D_T = 1.25 * ureg.inch
r_Cl = 2.65E6 * ureg.mg / ureg.L
r_Co = 1.138E6 * ureg.mg / ureg.L
P_HD = 0.1

def C_Prec(Dose: 'in mg/L'): # Function for calculating precipitate concentration
	return Dose.to(ureg.mg/ureg.L)*1.3/0.027/13 # Changed 1.039384 to 1.3

def phi_0(Dose: 'in mg/L',Inf: 'in NTU'): # Function for calculation phi_0
	x = C_Prec(Dose)/r_Co + Inf.to(ureg.mg/ureg.L)/r_Cl
	return x.to(ureg.dimensionless) 

def P_ClSphere(P_HD): # Function for SA to volume ratio for clay normalized to that of a sphere
	return (0.5 + P_HD)*(2/(3*P_HD))**(2/3)

def P_AClATot(Inf: 'in NTU',D_Cl: 'in um',D_T: 'in inches',P_HD): #Function to normalize surface area of clay by total surface area (clay + walls)
	x = 1/(1+(2*D_Cl.to(ureg.m)/(3*D_T.to(ureg.m)*P_ClSphere(P_HD)*Inf/r_Cl)))
	return x.to(ureg.dimensionless)

def Gamma(Inf,D,D_T,D_Cl,P_HD,D_Co): # Function for calculating Gamma
	return 1 - np.exp((-phi_0(D,0*ureg.mg/ureg.L)*D_Cl)/(phi_0(0*ureg.mg/ureg.L,Inf)*D_Co)*1/np.pi*(P_AClATot(Inf,D_Cl,D_T,P_HD))/(P_ClSphere(P_HD)))

def D_Sep(phi,D_Cl): # Function to give the particle separation distance.
	x = (np.pi/6/phi)**(1/3)*D_Cl
	return x.to(ureg.micrometer)

def P_V(Gam,t,EDR,nu,phi):
	x = (2*Gam-Gam**2)*t*np.sqrt(EDR.to(ureg.m**2/ureg.s**3)/nu.to(ureg.m**2/ureg.s))*phi**(2/3)
	return x.to(ureg.dimensionless)

def P_I(Gam,t,EDR,D_Cl,phi):
	x = (2*Gam-Gam**2)*t*(EDR.to(ureg.m**2/ureg.s**3)/D_Cl.to(ureg.m)**2)**(1/3)*phi**(8/9)
	return x.to(ureg.dimensionless)

def pC_I(N, k):
	return 9.0/8.0*np.log10(8.0/9.0*np.pi*k*N*(6.0/np.pi)**(8.0/9.0) + 1)

def pC_V(N, k): 
	return 1.5*np.log10(2.0/3.0*np.pi*k*N*(6.0/np.pi)**(2.0/3.0) + 1)

def Lambda_V(t,D,Inf,k_V,Gam,EDR,nu,phi,D_Cl):
	x = D_Sep(phi_0(D,Inf*10**(-pC_V(P_V(Gam,t,EDR,nu,phi),k_V))),D_Cl)
	return x.to(ureg.micrometer)

def Lambda_I(t,D,Inf,k_I,Gam,EDR,phi,D_Cl):
	x = D_Sep(phi_0(D,Inf*10**(-pC_I(P_I(Gam,t,EDR,D_Cl,phi),k_I))),D_Cl)
	return x.to(ureg.micrometer)

# Calculate New Gammas

Gam_7_90 = Gamma(Inf,D,D_T,D_Cl,P_HD,D_Co)

PV_7_90 = PV / Gam * (2*Gam_7_90-Gam_7_90**2)
PI_7_90 = PI / Gam * (2*Gam_7_90-Gam_7_90**2)

Gam_7_20 = Gamma(Inf,D,D_T,D_Cl,P_HD,20*ureg.nm)

PV_7_20 = PV / Gam * Gam_7_20
PI_7_20 = PI / Gam * Gam_7_20

# Fit Curves

N_c = np.arange(0.5,1.5E3,0.5)

k_I, k_Ivar = curve_fit(pC_I, PI_7_90, pC)
k_V, k_Vvar = curve_fit(pC_V, PV_7_90, pC)

# Make Figures
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif',size='14')
plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]
params = {'text.usetex' : True,
          'font.size' : 14,
          'font.family' : 'serif',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

plt.figure(0)
plt.semilogx(PI_7_90[0:16],pC[0:16],'ro',ms=10,label=r'Data')
plt.semilogx(PI_7_90[16:19],pC[16:19],'ro',ms=10,markerfacecolor='none')
plt.semilogx(N_c,pC_I(N_c, k_I),'r',ms=10,label=r'Inertial Model')
plt.xlabel(r'$\overline{\alpha}\theta\left(\frac{\overline{\varepsilon}}{\overline{d}_{\mathrm{P}}^2}\right)^{1/3}\phi_0^{8/9}$')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0.5,1.5E3,0,3.0])
#plt.title('Inertial Fit')
plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Inertial_Paper.png',format='png')
plt.savefig('Inertial_Paper.eps',format='eps')

plt.figure(1)
plt.semilogx(PV_7_90[0:16],pC[0:16],'bo',ms=10,label=r'Data') 
plt.semilogx(PV_7_90[16:19],pC[16:19],'bo',ms=10,markerfacecolor='none') 
plt.semilogx(N_c,pC_V(N_c,k_V),'b',ms=10,label=r'Viscous Model') 
plt.xlabel(r'$\overline{\alpha}\theta \overline{G}\phi_0^{2/3}$')
plt.ylabel(r'$\mathrm{p}C^{*}$')
plt.axis([0.5,1.5E3,0,3.0])
#plt.title('Viscous Fit')
plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Viscous_Paper.png',format='png')
plt.savefig('Viscous_Paper.eps',format='eps')
'''
D_Plot = np.arange(0.01,200,0.01)*ureg.mg/ureg.L

GamPlot_7_90 = Gamma(np.mean(Inf),D_Plot,D_T,D_Cl,P_HD,D_Co)
GamPlot_7_20 = Gamma(np.mean(Inf),D_Plot,D_T,D_Cl,P_HD,20*ureg.nm)

plt.figure(2)
plt.semilogx(D_Plot,GamPlot_7_20,'r',ms=10,label=r'$D_\mathrm{PACl}=20\:\mathrm{nm}$, $D_\mathrm{Clay}=7\:\mathrm{\muup m}$') 
plt.semilogx(D_Plot,GamPlot_7_90,'b',ms=10,label=r'$D_\mathrm{PACl}=90\:\mathrm{nm}$, $D_\mathrm{Clay}=7\:\mathrm{\muup m}$') 
plt.xlabel(r'Coagulant Dose (mg/L as Al)')
plt.ylabel(r'$\overline{\Gamma}$')
plt.axis([0.01,2E2,0,1.0])
#plt.title('Comparison of Gammas')
plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Gammas_New.png',format='png')
plt.savefig('Gammas_New.eps',format='eps')

t_Plot = np.arange(0,t[13]/ureg.s,1)*ureg.s
Kolm13 = (nu[13]**3/EDR[13])**(1/4)
K_Plot = np.ones(len(t_Plot))*Kolm13.to(ureg.micrometer)

plt.figure(4)
plt.semilogx(t_Plot,Lambda_V(t_Plot,D[13],Inf[13],k_V,Gam_7_90[13],EDR[13],nu[13],phi[13],D_Cl),'b',ms=10,label=r'Viscous')
plt.semilogx(t_Plot,Lambda_I(t_Plot,D[13],Inf[13],k_I,Gam_7_90[13],EDR[13],phi[13],D_Cl),'r',ms=10,label=r'Inertial')
plt.semilogx(t_Plot,K_Plot,'k--',ms=10,label=r'Kolmogorov Microscale')
plt.xlabel(r'Time in Flocculator (s)')
plt.ylabel('$\overline{\Lambda} \: (\muup \mathrm{m})$')
plt.axis([1,t[13]/ureg.s,60,180])
#plt.title('Comparison of Separation Distances')
plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Lambdas_New.png',format='png')
plt.savefig('Lambdas_New.eps',format='eps')

plt.figure(5)
plt.plot(D,Eff,'k+',ms=10)
plt.xlabel(r'PACl Dose (mg/L as Al)')
plt.ylabel(r'Effluent Turbidity (NTU)')
#plt.axis([1,t[13]/ureg.s,60,180])
#plt.title('Comparison of Separation Distances')
#plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Effluent.png',format='png')
plt.savefig('Effluent.eps',format='eps')
'''

# Calculate r^2 values.

#r2_I = r2_score(pC[0:14], pC_I(PI_7_90[0:14],k_I))
#r2_V = r2_score(pC[0:14], pC_V(PV_7_90[0:14],k_V))

# Calculate SSE & RMSE.

def SSE(A,B):
	return np.sum((A-B)**2)

def RMSE(A,B):
	return np.sqrt(np.mean((A-B)**2))

# Whole Dataset
SSE_I = SSE(pC, pC_I(PI_7_90,k_I))
SSE_V = SSE(pC, pC_V(PV_7_90,k_V))

RMSE_I = RMSE(pC, pC_I(PI_7_90,k_I))
RMSE_V = RMSE(pC, pC_V(PV_7_90,k_V))

# Last two points removed.
k_Is, k_Ivars = curve_fit(pC_I, PI_7_90[0:16], pC[0:16])
k_Vs, k_Vvars = curve_fit(pC_V, PV_7_90[0:16], pC[0:16])

SSE_Is = SSE(pC[0:16], pC_I(PI_7_90[0:16],k_Is))
SSE_Vs = SSE(pC[0:16], pC_V(PV_7_90[0:16],k_Vs))

RMSE_Is = RMSE(pC[0:16], pC_I(PI_7_90[0:16],k_Is))
RMSE_Vs = RMSE(pC[0:16], pC_V(PV_7_90[0:16],k_Vs))

print('\nInertial Values\n')

print('Sum Squared Error: ',str(SSE_I))
print('Root Mean Squared Error: ',str(RMSE_I))
print('k: ',str(k_I))
print('Subset Sum Squared Error: ',str(SSE_Is))
print('Subset Root Mean Squared Error: ',str(RMSE_Is))
print('Subset k: ',str(k_Is))

print('\nViscous Values\n')

print('Sum Squared Error: ',str(SSE_V))
print('Root Mean Squared Error: ',str(RMSE_V))
print('k: ',str(k_V))
print('Subset Sum Squared Error: ',str(SSE_Vs))
print('Subset Root Mean Squared Error: ',str(RMSE_Vs))
print('Subset k: ',str(k_Vs))

plt.close('all')
