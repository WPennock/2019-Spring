# Individual Data Processor
This script will process single experiments and export processed data to a .csv file which can then be read by the aggregate data processor.

## Imports
```python
from aguaclara.play import *
import aguaclara.research.floc_model as floc
import aguaclara.research.procoda_parser as pro
import scipy.stats as stats
import doctest
import pdb
import pandas as pd

plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]
params = {'text.usetex' : True,
          'font.size' : 10,
          'font.family' : 'serif',
          'text.latex.unicode': True,
          'axes.facecolor': 'white',
          'axes.labelcolor': 'black',
          'savefig.facecolor': 'white',
          'axes.edgecolor': 'black',
          'savefig.edgecolor': 'black',
          'xtick.color': 'black',
          'ytick.color': 'black',
          'grid.linestyle': 'None'
          }
plt.rcParams.update(params)
```

## Load Data
```python
MetaID = 88
# Drive depends on computer where data are being processed
Drive = "C:\\Users\\whp28\\Google Drive\\AGUACLARA DRIVE\\AguaClara Grads\\William Pennock\\"
# Drive = "C:\\Users\\William\\Google Drive\\AGUACLARA DRIVE\\AguaClara Grads\\William Pennock\\"
Meta_Path = Drive + "Meta File.xls"
Meta = pd.read_excel(Meta_Path)
Path = Meta.loc[MetaID-1, "C:\\Users\\whp28\\Google Drive\\AGUACLARA DRIVE\\AguaClara Grads\\William Pennock\\"] + "\\"
Day = Meta.loc[MetaID-1, "Begin"]
state_file = Drive + Path + "statelog " + Day + ".xls"
data_file = Drive + Path + "datalog " + Day + ".xls"
StateFile = pd.read_csv(state_file,'\t')
DayFrac = pro.column_of_data(data_file, 0, 0)
StateFile
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]
StateFile["Day fraction since midnight on " + Day.replace("-","/")]
Start = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 1].index].values
StartLoc = np.where(DayFrac==find_nearest(DayFrac,Start))[0][0]

End = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 0].index].values
EndLoc = np.where(DayFrac==find_nearest(DayFrac,End))[0][0]
```

## Data Inspection
```python
# Process Variables
# Load the entirety of these variables, because they give a picture of the overall process.
Time = pro.column_of_time(data_file, 0)
State = pro.column_of_data(data_file, 0, 17)

Clay = pro.column_of_data(data_file, 0, 8)
SWaT = pro.column_of_data(data_file, 0, 9)

plt.clf(), plt.close('all')
fig0 = plt.figure(0)
ax1 = fig0.add_subplot(111)
line1 = ax1.plot(Time[StartLoc:EndLoc],State[StartLoc:EndLoc],'r',label="State")
plt.ylabel("State")
plt.xlabel("Time (day)")
ax2 = fig0.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time,Clay,'g',label="Clay Pump Fraction")
line3 = ax2.plot(Time,SWaT,'b',label="SWaT Pump Fraction")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Pump Fraction")
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=0)
plt.savefig("State and Pump Fraction "+str(MetaID)+".png")
plt.show()

pH = pro.column_of_data(data_file, 0, 1)

T_Air = pro.column_of_data(data_file, 0, 3)*u.degC
T_Wat = pro.column_of_data(data_file, 0, 4)*u.degC
T_AC = pro.column_of_data(data_file, 0, 5)*u.degC

plt.clf(), plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
line1 = ax1.plot(Time,pH,'r',label="pH")
plt.xlabel("Time (day)")
plt.ylabel("pH")
ax2 = fig1.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time[StartLoc:EndLoc],T_Air[StartLoc:EndLoc],'g',label="Air Temperature")
line3 = ax2.plot(Time[StartLoc:EndLoc],T_Wat[StartLoc:EndLoc],'xkcd:purple',label="Water Temperature")
line4 = ax2.plot(Time[StartLoc:EndLoc],T_AC[StartLoc:EndLoc],'b',label="AC Temperature")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Temperature (C)")
plt.ylim([20,26])
# ax2.legend()
lines = line1 + line2 + line3 + line4
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=3)
plt.savefig("pH and Temperature "+str(MetaID)+".png")
plt.show()

hL = pro.column_of_data(data_file, 0, 2)*u.cm
d_AC = pro.column_of_data(data_file, 0, 6)*u.cm

Bal = pro.column_of_data(data_file, 0, 7)*u.g

plt.clf(), plt.close('all')
fig2 = plt.figure(2)
ax1 = fig2.add_subplot(111)
line1 = ax1.plot(Time[StartLoc:EndLoc],Bal[StartLoc:EndLoc],'g',label="Balance")
plt.xlabel("Time (day)")
plt.ylabel("Mass (g)")
ax2 = fig2.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time[StartLoc:EndLoc],hL[StartLoc:EndLoc],'r',label="Head Loss")
line3 = ax2.plot(Time[StartLoc:EndLoc],d_AC[StartLoc:EndLoc],'b',label="AC Height")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Height (cm)")
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=0)
plt.savefig("Balance and Water Levels "+str(MetaID)+".png")
plt.show()

# Process Data

Inf = pro.column_of_data(data_file, 0, 11)*u.NTU
Eff = pro.column_of_data(data_file, 0, 12)*u.NTU
Abs = pro.column_of_data(data_file, 0, 19)

plt.clf(), plt.close('all')
fig3 = plt.figure(3)
ax1 = fig3.add_subplot(111)
line1 = ax1.plot(Time[StartLoc:EndLoc],Abs[StartLoc:EndLoc],'g',label="UV 254")
plt.xlabel("Time (day)")
plt.ylabel("Absorbance")
plt.ylim([-2,2])
ax2 = fig3.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time[StartLoc:EndLoc],Inf[StartLoc:EndLoc],'r',label="Influent")
line3 = ax2.plot(Time[StartLoc:EndLoc],Eff[StartLoc:EndLoc],'b',label="Effluent")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Turbidity (NTU)")
plt.ylim([0,200])
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=0)
plt.savefig("Absorbance and Turbidity "+str(MetaID)+".png")
plt.show()
```


## Dividing Experiments
```python
# Finding indices of experiments

Begins = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 2].index].values
BeginData = [find_nearest(DayFrac,l) for l in Begins]
BeginLoc = [np.where(DayFrac==l)[0][0] for l in BeginData]
BeginLoc
BeginLoc.append(EndLoc)
BeginLoc

Reads = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 3].index].values
ReadData = [find_nearest(DayFrac,l) for l in Reads]
ReadLoc = [np.where(DayFrac==l)[0][0] for l in ReadData]
ReadLoc

Flush = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 4].index].values
FlushLoc = np.where(DayFrac==find_nearest(DayFrac,Flush))[0][0]
ReadLoc.append(FlushLoc)

# Calculate lag time
def Q_c(v_c, D, L, alpha):
    return (np.pi/4*D**2*v_c*(L/D*np.cos(alpha)+np.sin(alpha))).to(u.mL/u.min)
L_SWaT = 86*u.cm
D_SWaT = 1.049*u.inch
a_SWaT = 60*u.deg    
v_c = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])*u.mm/u.s
Q_SWaT = Q_c(v_c,D_SWaT,L_SWaT,a_SWaT)
Q_SWaT
V_SWaT = np.pi*D_SWaT**2/4*L_SWaT
T_SWaT = (V_SWaT/Q_SWaT).to(u.s)
T_SWaT
V_Plant = (56.35*u.m*np.pi/4*(1.25*u.inch)**2).to(u.L)
T_Plant = Meta.loc[MetaID-1, "Residence Time (s)"]*u.s
T_Total = T_Plant + T_SWaT
Lag = [int(np.round(l/(5*u.s))) for l in T_Total]

# Inspect performance of different runs
plt.clf(),plt.close('all')
fig4 = plt.figure(4)
for i in range(0,len(ReadLoc)):
    ax = fig4.add_subplot(2,3,i+1)
    ax.plot(Time[ReadLoc[i]:BeginLoc[i+1]],Inf[ReadLoc[i]-Lag[i]:BeginLoc[i+1]-Lag[i]],'r-')
    ax.plot(Time[ReadLoc[i]:BeginLoc[i+1]],Eff[ReadLoc[i]:BeginLoc[i+1]],'b-')
plt.savefig("Experimental Performance "+str(MetaID)+".png")
plt.show()
# Adjust Endpoints
## Subplot 2
plt.clf(),plt.close('all')
plt.figure(5)
plt.plot(Time[ReadLoc[1]:BeginLoc[1+1]],Inf[ReadLoc[1]-Lag[1]:BeginLoc[1+1]-Lag[1]],'r-')
plt.plot(Time[ReadLoc[1]:BeginLoc[1+1]],Eff[ReadLoc[1]:BeginLoc[1+1]],'b-')
# plt.axis([0.057,0.061,70,75])
plt.show()

ReadLoc[1] = np.where(Time.magnitude == find_nearest(Time.magnitude,0.0585))[0][0]

## Subplot 3
plt.clf(),plt.close('all')
plt.figure(5)
plt.plot(Time[ReadLoc[2]:BeginLoc[2+1]],Inf[ReadLoc[2]-Lag[1]:BeginLoc[2+1]-Lag[1]],'r-')
plt.plot(Time[ReadLoc[2]:BeginLoc[2+1]],Eff[ReadLoc[2]:BeginLoc[2+1]],'b-')
# plt.axis([0.057,0.061,70,75])
plt.show()

## Subplot 6
plt.clf(),plt.close('all')
plt.figure(6)
plt.plot(Time[ReadLoc[5]:BeginLoc[5+1]],Inf[ReadLoc[5]:BeginLoc[5+1]],'r-')
plt.plot(Time[ReadLoc[5]:BeginLoc[5+1]],Eff[ReadLoc[5]:BeginLoc[5+1]],'b-')
# plt.axis([0.0952,0.0967,81,84])
plt.show()

T_SWaT[5].to(u.day)
ReadLoc[5] = np.where(Time.magnitude == find_nearest(Time.magnitude,0.0952))[0][0]
BeginLoc[5+1] = np.where(Time.magnitude == find_nearest(Time.magnitude,0.0967))[0][0]
```



## Find Average Values
```python
# Lagged Values
pH_avg = np.zeros(len(v_c))
Inf_avg = np.zeros(len(v_c))
hL_Floc = np.zeros(len(v_c))
T_AC_avg = np.zeros(len(v_c))
d_AC_avg = np.zeros(len(v_c))

for i in range(0,len(v_c)):
    pH_avg[i] = np.mean(pH[ReadLoc[i]-Lag[i]:BeginLoc[i+1]-Lag[i]])
    Inf_avg[i] = np.mean((Inf[ReadLoc[i]-Lag[i]:BeginLoc[i+1]-Lag[i]]).to(u.NTU)).magnitude
    hL_Floc[i] = np.mean((hL[ReadLoc[i]-Lag[i]:BeginLoc[i+1]-Lag[i]]).to(u.cm)).magnitude
    d_AC_avg[i] = np.mean((d_AC[ReadLoc[i]-Lag[i]:BeginLoc[i+1]-Lag[i]]).to(u.cm)).magnitude
  # Could adjust lag on T_AC to reflect time between tank and inlet  
    T_AC_avg[i] = np.mean((T_AC[ReadLoc[i]-Lag[i]:BeginLoc[i+1]-Lag[i]]).to(u.degC)).magnitude

# Real Time Values
Abs_avg = np.zeros(len(v_c))
Eff_avg = np.zeros(len(v_c))
T_Wat_avg = np.zeros(len(v_c))
T_Air_avg = np.zeros(len(v_c))

for i in range(0,len(v_c)):
    Abs_avg[i] = np.mean(pH[ReadLoc[i]:BeginLoc[i+1]])
    Eff_avg[i] = np.mean((Eff[ReadLoc[i]:BeginLoc[i+1]]).to(u.NTU)).magnitude
    T_Wat_avg = np.mean((T_Wat[ReadLoc[i]:BeginLoc[i+1]]).to(u.degC)).magnitude
    T_Air_avg = np.mean((T_Air[ReadLoc[i]:BeginLoc[i+1]]).to(u.degC)).magnitude

# Constant values
Q = [((V_Plant/T_Plant).to(u.mL/u.s)).magnitude for l in range(0,len(v_c))]*u.mL/u.s # Chemical addition all <1% of total, SWaT varies from 1-4%.
v_c = (v_c.to(u.mm/u.s)).magnitude
ID = [MetaID + l for l in range(0,len(v_c))] # Chemical addition all <1% of total, SWaT varies from 1-4%.

plt.clf(),plt.close('all')
Time_0 = -25
plt.figure(7)
plt.plot(Time[Time_0:len(Time)-6],Abs[Time_0:len(Time)-6],'g')
plt.plot(Time[Time_0:len(Time)-6],pH[Time_0:len(Time)-6],'r')
plt.plot(Time[Time_0:len(Time)-6],Inf[Time_0:len(Time)-6],'b')
plt.show()

Abs_0 = [np.mean(Abs[Time_0:len(Time)-6]) for l in range(0,len(ReadLoc))]
pH_0 = [np.mean(pH[Time_0:len(Time)-6]) for l in range(0,len(ReadLoc))]
Inf_0 = [np.mean(Inf[Time_0:len(Time)-6]) for l in range(0,len(ReadLoc))]
hL_0 =
plt.plot(Time, hL)
```
## Determining Coagulant Dose
```python
C_PS = Meta.loc[MetaID-1, "Coagulant\nStock Conc. (mg/L)"]*u.mg/u.L
rho_PSS = 1.27*u.kg/u.L # Assumed from prior experimental data
rho_H2O = pc.density_water(24*u.degC)
rho_H2O
rho_PSS.to(u.kg/u.m**3)
C_PSS = 71.12*u.g/u.L
Ratio_P = (C_PS/C_PSS).to(u.dimensionless)
rho_P = Ratio_P*rho_PSS + (1-Ratio_P)*rho_H2O
rho_P

plt.clf(),plt.close('all')
fig4 = plt.figure(4)
for i in range(0,len(ReadLoc)):
    ax = fig4.add_subplot(2,3,i+1)
    ax.plot(Time[ReadLoc[i]:BeginLoc[i+1]],Bal[ReadLoc[i]-Lag[i]:BeginLoc[i+1]-Lag[i]],'r-')
# plt.savefig("Balance "+str(MetaID)+".png")
plt.show()

linreg = np.zeros([len(ReadLoc),5])
for i in range(0,len(ReadLoc)):
    linreg[i][:] = stats.linregress(Time,Bal)
slope_P = linreg[:,0]
int_P = linreg[:,1]
r_value_P = linreg[:,2]
r_value_P**2    

Dose = (-slope_P*u.g/u.day)/(rho_P)/Q*C_PS
Dose.to(u.mg/u.L)
Dose = np.zeros(len(ReadLoc))
for i in range(0,len(ReadLoc)):
    Dose[i] = (slope_P*u.g/u.day)*
Dose  
```

## Calculating Head Loss
```python
hL_Meas = Meta.loc[MetaID-1, "Head Loss (cm)"]*u.cm
plt.plot(Time[EndLoc+30:-20],d_AC[EndLoc+30:-20])
Datum = np.round(np.mean(d_AC[EndLoc+30:-20]),2)
Datum
hL_Total = [np.mean(hL_Meas + d_AC[ReadLoc[i]:BeginLoc[i+1]] - Datum) for l in range(0,len(ReadLoc))]
```


## Exporting Data
```python
Matrix = np.column_stack()
```
