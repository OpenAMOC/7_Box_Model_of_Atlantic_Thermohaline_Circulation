# 7_Box_Model_of_Atlantic_Thermohaline_Circulation
An intermediate complexity model with 4 ocean boxes and 3 atmospheric boxes simulating the thermohaline circulation in the Atlantic.
For more info: https://www.openamoc.site/for-python

Once unzipped, you will find 3 python scripts. The core script you should run is AMOC_Main.py, with AMOC_Constants.py and AMOC_Bmodel.py in the same folder.
Under default settings the programe takes approx 5 and a half minutes to run (dependent on system capabilities). 
In the current folder two timeseries graphs are saved as .png files and an Excel document of the model's results. This will be overwritten with each run.

Variables:
tos,sos = Temp, Salinity of north box

tom,som = Temp, Salinity of surface current

ton,son = Temp, Salinity of south box

tod,sod= Temp, Salinity of deep current

tas, tam, tan= Temp for atmosphere over north Box, surface current and south box
