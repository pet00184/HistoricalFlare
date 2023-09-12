### GOES XRS Exploratory Data Analysis:

**Contents**:

`GOES_XRS_historical.fits`: FITS file that includes data for all GOES flares from 2017-August 2023. This is the .FITS file to be used when performing the parameter search. Depending on size, either new columns will be added when more arrays are made, or a new FITS with just those arrays (e.g. derivatives) will be incorporated.

`FITS_plots`: visualization of all the flares.

`parameter_search.py`: source code for the parameter search.

`plotting_results.py`: source code for the plotting of results done after the parameter search.

`run_paramsearch.py`: code to be run for the parameter search! As new parameters are made, new functions will be added to this file.

**To Run:**
To run an example search (using the flux value for XRSB), run `python3 run_paramsearch.py` in this folder
