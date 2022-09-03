This file describes how to obtain the datasets from ECB SDW and Eurostat

*********************************************************
ECB SDW
***Obtain prices indices

- https://sdw.ecb.europa.eu/

Type these reference number in the search bar:
"ICP.M.FR.N.000000.4.INX", "ICP.M.FR.N.010000.4.INX", "ICP.M.FR.N.020000.4.INX", "ICP.M.FR.N.030000.4.INX", "ICP.M.FR.N.040000.4.INX", "ICP.M.FR.N.050000.4.INX", "ICP.M.FR.N.060000.4.INX", "ICP.M.FR.N.070000.4.INX", "ICP.M.FR.N.080000.4.INX", "ICP.M.FR.N.090000.4.INX", "ICP.M.FR.N.100000.4.INX", "ICP.M.FR.N.110000.4.INX", "ICP.M.FR.N.120000.4.INX"

These are the references for wide categories price indexes and correspond to:
'HICP - Overall index',	'HICP - FOOD','HICP - ALCOHOL, TOBACCO','HICP - CLOTHING','HICP - ELECTRICITY, GAS','HICP - FURNISHINGS','HICP - HEALTH','HICP - TRANSPORT','HICP - COMMUNICATION','HICP - RECREATION ','HICP - EDUCATION ','HICP - RESTAURANTS','HICP - MISCELLANEOUS'

*********************************************************
EUROSTAT
***Obtain Final consumption expenditure share by purpose

- https://ec.europa.eu/eurostat/databrowser/view/NAMA_10_CO3_P3__custom_3301056/default/table?lang=en

Final consumption expenditure of households by consumption purpose (COICOP 3 digit)
online data code: NAMA_10_CO3_P3
Source of data: Eurostat

- Select Euro Area 19 from 2015 and Euro Area before.
- Select the 12 purposes of consumption and indicate that you would like to show shares.

*********************************************************
EUROSTAT
***Obtain consumption expenditure structure by quintile:
https://ec.europa.eu/eurostat/databrowser/view/HBS_STR_T223__custom_3301043/default/table?lang=en

Structure of consumption expenditure by income quintile and COICOP consumption purpose (online data code: HBS_STR_T223 )
Source of data: Eurostat

- Select Euro Area 19 from 2015 and Euro Area before.
- Select the 12 purposes of consumption and also indicate the 5 quintiles.

- Values obtained from Eurostat are in value. To obtain shares, I sum per year the total consumption expenditure across quintiles and divide initial value per purpose by this total consumption expenditure.