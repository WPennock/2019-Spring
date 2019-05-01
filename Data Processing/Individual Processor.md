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
MetaID = 40
# Drive depends on computer where data are being processed
# Drive = "C:\\Users\\whp28\\Google Drive\\AGUACLARA DRIVE\\AguaClara Grads\\William Pennock\\
Drive = "C:\\Users\\William\\Google Drive\\AGUACLARA DRIVE\\AguaClara Grads\\William Pennock\\"
Meta_Path = Drive + "Meta File.xls"
Meta = pd.read_excel(Meta_Path)
Path = Meta.loc[MetaID-1, "C:\\Users\\whp28\\Google Drive\\AGUACLARA DRIVE\\AguaClara Grads\\William Pennock\\"] + "\\"
Day = Meta.loc[MetaID-1, "Begin"]
state_file = Drive + Path + "statelog " + Day + ".xls"
data_file = Drive + Path + "datalog " + Day + ".xls"
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
line1 = ax1.plot(Time,State,'r',label="State")
plt.ylabel("State")
ax2 = fig0.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time,Clay,'g',label="Clay Pump Fraction")
line3 = ax2.plot(Time,SWaT,'b',label="SWaT Pump Fraction")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Pump Fraction")
# ax2.legend()
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=0)
plt.show()

pH = pro.column_of_data(data_file, 0, 1)

T_Air = pro.column_of_data(data_file, 0, 3)*u.degC
T_Wat = pro.column_of_data(data_file, 0, 4)*u.degC
T_AC = pro.column_of_data(data_file, 0, 5)*u.degC

plt.clf(), plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
line1 = ax1.plot(Time,pH,'r',label="pH")
plt.ylabel("pH")
ax2 = fig1.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time,T_Air,'g',label="Air Temperature")
line3 = ax2.plot(Time,T_Wat,'tab:purple',label="Water Temperature")
line4 = ax2.plot(Time,T_AC,'b',label="AC Temperature")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Temperature (C)")
plt.ylim([20,30])
# ax2.legend()
lines = line1 + line2 + line3 + line4
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels)
plt.show()

hL = pro.column_of_data(data_file, 0, 2)*u.cm
d_AC = pro.column_of_data(data_file, 0, 6)*u.cm

Bal = pro.column_of_data(data_file, 0, 7)*u.g

# Process Data

Inf_All = pro.column_of_data(data_file, 0, 11)
Eff_All = pro.column_of_data(data_file, 0, 12)
Abs_All = pro.column_of_data(data_file, 0, 19)

plt.figure(0)
plt.plot(Time, Eff_All)
plt.ylim([0,100])
```
## Averaging Data
```python
# Performance Variables
# Load these only in the reading state
Eff_Data = pro.read_state("4-27-2019", 3, 12, units='NTU', path=Path)
Eff_Data
plt.plot(Eff_Data[0],Eff_Data[1])
plt.show()
v_c =
Inf = pro.average_state("4-27-2019",3,11,"NTU",Path)
Eff = pro.average_state("4-27-2019",3,12,"NTU",Path)
Abs = pro.average_state("4-27-2019",3,19,"NTU",Path)
Eff
# Performance Variable Adjustment
# For cases where the data need to be adjusted.

```

## Determining Coagulant Dose
```python

```

## Exporting Data
```python
```
