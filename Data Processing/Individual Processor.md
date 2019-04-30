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
```

## Load Data
```python
MetaID = 40
Meta_Path = r"C:\\Users\William\Google Drive\AGUACLARA DRIVE\AguaClara Grads\William Pennock\Meta File.xls"
Meta = pd.read_excel(Meta_Path)
Meta.rename(index={0,1})
Path = Meta.loc[MetaID-1, "C:\\Users\\whp28\\Google Drive\\AGUACLARA DRIVE\\AguaClara Grads\\William Pennock\\"]
Datalog = 
```
## Bulk Data
```python
# Process Variables
# Load the entirety of these variables, because they give a picture of the overall process.
Time = pro.column_of_time(Path, 0)
pH = pro.column_of_data(Path, 1, 1)
hL =
```
##
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

## Determine Coagulant Dose
```python

```
