# ELEvoHI
ELlipse Evolution model based on Heliospheric Imager observations. ELEvoHI predicts arrival times and speeds of coronal mass ejections at any target in the inner heliosphere. Detailed information on the model can be found in T. Rollett et al. (2016, single-run version) and T. Amerstorfer et al. (2018, ensemble version).

In the following, the three ambient solar wind models currently implemented in ELEvoHI are described.


**HUX:**

HUX provides a static solution of the solar wind speed for a full Carrington rotation in the solar equatorial plane.
The data is given in a simple 'txt' format while the initial date and time is given in an additional file.
Data description:
    Radial extent:  5 - 430 RSun (1 RSun resolution)
    Longitude:      2° (with the longituded of Earth set to 60°)
    
References: Reiss, M. et al. 2019 ApJS, 240, 35; Reiss, M. et al. 2020 ApJ, 891, 165


**HUXt:**

HUXt provides a time dependent solution for the ambient solar wind in the ecliptic plane. The data is provided in a 'hdf5' format. 
The initial start time for this model is set in an individual 'txt' file.
The structure of the file has to contain the following:
    'lon':          longitudinal extent
    'dlon':         resolution of the longitude
    'r':            radial extent
    'cr_lon_init':  initial Carrinton rotation longitude
    'dt_out':       temporal resolution
    'v_grid':       solar wind speed array
    'time_out':     time array

Data description:
    Radial extent:  21.5 - 300.5 RSun (1 RSun resolution)
    Longitude:      0.7° resolution
    Time:           3.865 seconcds resolution 


**EUHFORIA**

EUHFORIA provides a static solution of the solar wind speed and density for a full Carrington rotation and for different latitudes.
The data is saved in 'h5' format and the initial start time of the model output is saved in an individual 'txt' file.
The structure of the file has to contain the following:
    'lon_edges':    longitude (of the edges of the data grid)
    'r_edges':      radius (of the edges of the data grid)
    'vr_lat_m6':    ambient solar wind speed for latitude 6°
    'n_lat_m6':     ambient solar wind density for latitude 6°

Data description:
    Radial extent:  20.56 - 324.43 RSun (0.94 RSun resolution)
    Longitude:      1° resolution
    
References: Pomoell, J. & Poedts, S., 2018 SWSC, 8, A35
