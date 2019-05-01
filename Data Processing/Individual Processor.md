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
StateFile = pd.read_csv(state_file,'\t')
DayFrac = pro.column_of_data(data_file, 0, 0)

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

Start = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 1].index].values
StartLoc = np.where(DayFrac==find_nearest(DayFrac,Start))[0][0]

End = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 0].index].values
EndLoc = np.where(DayFrac==find_nearest(DayFrac,End))[0][0]
```

## Data Inspection
```python
# Process Variables
# Load the entirety of these variables, because they give a picture of the overall process.
Time = pro.column_of_time(data_file, 0)[StartLoc:EndLoc]
State = pro.column_of_data(data_file, 0, 17)[StartLoc:EndLoc]

Clay = pro.column_of_data(data_file, 0, 8)[StartLoc:EndLoc]
SWaT = pro.column_of_data(data_file, 0, 9)[StartLoc:EndLoc]

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
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=0)
plt.show()

pH = pro.column_of_data(data_file, 0, 1)[StartLoc:EndLoc]

T_Air = pro.column_of_data(data_file, 0, 3)[StartLoc:EndLoc]*u.degC
T_Wat = pro.column_of_data(data_file, 0, 4)[StartLoc:EndLoc]*u.degC
T_AC = pro.column_of_data(data_file, 0, 5)[StartLoc:EndLoc]*u.degC

plt.clf(), plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
line1 = ax1.plot(Time,pH,'r',label="pH")
plt.ylabel("pH")
ax2 = fig1.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time,T_Air,'g',label="Air Temperature")
line3 = ax2.plot(Time,T_Wat,'xkcd:purple',label="Water Temperature")
line4 = ax2.plot(Time,T_AC,'b',label="AC Temperature")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Temperature (C)")
plt.ylim([20,26])
# ax2.legend()
lines = line1 + line2 + line3 + line4
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=1)
plt.show()

hL = pro.column_of_data(data_file, 0, 2)[StartLoc:EndLoc]*u.cm
d_AC = pro.column_of_data(data_file, 0, 6)[StartLoc:EndLoc]*u.cm

Bal = pro.column_of_data(data_file, 0, 7)[StartLoc:EndLoc]*u.g

plt.clf(), plt.close('all')
fig2 = plt.figure(2)
ax1 = fig2.add_subplot(111)
line1 = ax1.plot(Time,Bal,'g',label="Balance")
plt.ylabel("Mass (g)")
ax2 = fig2.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time,hL,'r',label="Head Loss")
line3 = ax2.plot(Time,d_AC,'b',label="AC Height")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Height (cm)")
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=0)
plt.show()
# Process Data

Inf_All = pro.column_of_data(data_file, 0, 11)
Eff_All = pro.column_of_data(data_file, 0, 12)
Abs_All = pro.column_of_data(data_file, 0, 19)

plt.clf(), plt.close('all')
fig3 = plt.figure(3)
ax1 = fig3.add_subplot(111)
line1 = ax1.plot(Time,Abs_All,'g',label="UV 254")
plt.ylabel("Absorbance")
plt.ylim([-2,2])
ax2 = fig3.add_subplot(111,sharex=ax1,frameon=False)
line2 = ax2.plot(Time,Inf_All,'r',label="Influent")
line3 = ax2.plot(Time,Eff_All,'b',label="Effluent")
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Turbidity (NTU)")
plt.ylim([0,200])
lines = line1 + line2 + line3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc=0)
plt.show()
```

## Dividing Experiments
```python
# Finding indices of experiments


Begins = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 2].index].values
BeginData = [find_nearest(DayFrac,l) for l in Begins]
BeginLoc = [np.where(DayFrac==l) for l in BeginData]



Reads = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 3].index].values
ReadData = [find_nearest(DayFrac,l) for l in Reads]
ReadLoc = [np.where(DayFrac==l) for l in ReadData]

Flush = StateFile["Day fraction since midnight on " + Day.replace("-","/")].iloc[ StateFile[StateFile[" State ID"] == 4].index].values
FlushLoc = np.where(DayFrac==find_nearest(DayFrac,Flush))
ReadLoc.append(FlushLoc)

# Divide up experimental data
Inf = np.zeros(len(ReadLoc))
Eff = np.zeros(len(ReadLoc))
Abs = np.zeros(len(ReadLoc))

for i in range(0,len(ReadLoc))

```
## Experiments
### Experiment 1
```python
```
### Experiment 2
```python
```
### Experiment 3
```python
```
### Experiment 4
```python
```
## Determining Coagulant Dose
```python

```

## Exporting Data
```python
```
