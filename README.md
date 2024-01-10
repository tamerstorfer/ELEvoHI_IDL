
# IMPORTANT!
# This version of ELEvoHI is mainly written in IDL.
# It is not maintained anymore.
# Please use the current version coded in Python:
https://github.com/tamerstorfer/ELEvoHI

# ELEvoHI
ELlipse Evolution model based on Heliospheric Imager observations. ELEvoHI predicts arrival times and speeds of coronal mass ejections at any target in the inner heliosphere. Detailed information on the model can be found in T. Rollett et al. (2016, single-run version) and T. Amerstorfer et al. (2018, ensemble version).

In the following, it is described how ELEvoHI can be used based on data from heliospheric images alone, i.e. no coronagraph data are needed.

# How to use ELEvoHI

**1. Measure the CME front in heliospheric images**

The first and mandatory step to use ELEvoHI is to measure the propagating CME within heliospheric imager data. This can either be done in time-elongation maps or directly in the (running difference) images. A very convenient tool to perform these measurements is the SATPLOT tool developed at JPL. The resulting time-elongation needs to be stored as an IDL .sav-file and should contain an IDL-structure:

        IDL> help, track, /str
        ** Structure <2637478>, 4 tags, length=736, data length=736, refs=1:
           TRACK_DATE      STRING    Array[30]
           ELON            FLOAT     Array[30]
           ELON_STDD       FLOAT     Array[30]
           SC              STRING    'A'

Usually, we measure a CME front at least 5 times and build a mean track in order to get a good representative of the CME evolution. Note, that the standard deviation derived from these measurements is not used within ELEvoHI. Due to the ensemble character of the model, the differences between the ensemble members are much larger than any differences due to the measurement error of the time-elongation track.

**Important:** *Since ELEvoHI is used to predict arrival times and speeds at targets in the ecliptic, it is necessary to perform the measurements of the CME front within the ecliptic!*

**2. Derive the direction of motion**

The propagation direction is the only parameter that needs to be provided to run ELEvoHI. When using SATPLOT to measure the CME front, the direction can directly be calculated by performing a fit to the time-elongation profile. This procedure is built-in and yields fits based on the Fixed-Phi (FPF), the Harmonic Mean and the Self-Similar fitting approaches. From that, we also get a prediction based on the assumptions made by these methods (constant propagation speed and certain shape within the ecliptic). For ELEvoHI, we use the resulting propagation direction relative to our HI-observer (either STEREO-A or B) from the FPF method.

**3. Setting up the ELEvoHI source file**

The ELEvoHI_input.txt file is stored in the Code folder and should be edited directly without producing a folder for the predicted event. This is done automatically.

In the source file, we need to define the following parameters:

Is the CME observed from STEREO-A or STEREO-B?

        ELEvoHI source file

        ________________
        HI observer (A or B)
        A

When did the CME appear within HI1 for the first time?

        ________________
        first observed by HI1
        20081212

Do you use your own time-elongation measurements (e.g. from SATPLOT) or do you use measurements from HELCATS (www.helcats-fp7.eu)? In the latter case, change "user-defined" to "helcats".
When you run ELEvoHI for the first time for a specific event, you should not state the start and end point of the DBM fit in advance. The model provides the possibility to chose these points manually when they are not stated in the source file.
Since ELEvoHI is assuming a drag-based propagation, it is necessary that any early acceleration due to the Lorentz-force is excluded to be used for the prediction. Usually, a reasonable region to perform the DBM fit is from 30 to 100 solar radii. For some events, the initial distance can be also very close to the Sun or even farther out. This depends on the kinematics of each event individually. For chosing the start and end points manually (recommended for first run of each event), delete everything in the line after "user-defined".
In line 12, the path to your IDL .sav file containing the time-elongation track is given.

        ________________
        source of HI-tracks (helcats or user-defined/start/end of DBM fit)
        user-defined/5/20
        /home/tamerstorfer/hievents/SATPLOT/20081212_A/satplot_track.sav
        
Now, the parameters for the ensemble mode are set. The first is the inverse ellipse aspect ratio.
Usually, there is no need to change this parameter. We have realized that predicted flank hits often produce very late arrival times and suppose that this is due to the highly curved flanks when assuming a circular frontal shape. In order to extenuate this effect, it is recommended to vary the inverse aspect ratio between e.g. 0.6-0.9.
The first (0.7) and the second (1) numbers state the range in which the inverse aspect ration is varied. The third number (0.1) gives the step size.

        ________________
        inverse ellipse aspect ratio, f (b/a), ensemble: (f_min/f_max/d_f)
        0.7/1/0.1
     
In line 19, the propagation direction relative to the HI-observer (as derived from e.g. FPF) should be given. Note, that it cannot be smaller that 1!

        ________________
        direction of motion of CME apex, phi, ensemble: (phi_min/phi_max/d_phi)
        (angle between HI-observer and CME apex)
        46/66/2

For the angular half width within the ecliptic plane, a range between 30 and 50 degree and a step size of 5 degree seems to be a good choice.
The last value (kappa) defines the latitudinal extent of the CME. It is also defined as the half width of the CME. It is assumed not to change during propagation
Change this value according to your preferences.

        ________________
        angular half width, lambda, ensemble: (lambda_min/lambda_max/d_lambda/kappa)
        30/50/5/15

In line 25, the most probable in situ arrival target can be given. This is needed for the version of ELEvoHI, were in situ data from 1 AU is used as an approximation of the ambient solar wind speed that is influencing the CME evolution. It is not recommended to use this option and in the future, it will be deleted.

        ________________
        most likely in situ detection at 1 AU (Earth or A or B)
        Earth
        
The last change that can be made in the source file, is the actual in situ arrival target and time. When this information is provided, the difference between the predicted and detected arrival time is given in the results file and displayed in the figures. If this information is unknown (e.g. in case of a real-time prediction), leave these lines (28 and 29) blank.

        ________________
        in situ arrival s/c and times (MES, VEX, Earth, A, B)
        Earth
        2008-Dec-17 03:35
        
When the parameters needed are changed, save the file and start ELEvoHI as follows:

        elevohi, /save_results, /statistics, /silent, /forMovie, bgsw='stat'

The only interaction needed is the definition of the start and end points of the DBM fit. This is done by clicking into a figure and should be self-explanatory.

When ELEvoHI has finished, the result of the prediction is given in the terminal.

The Python-codes for the visualization of the prediction result can be found in https://github.com/tamerstorfer/ELEvoHI/tree/master/HI_animate_module.


To run ELEvoHI with a deformable front an ambient solar wind for the heliosphere is needed. Start ELEvoHI as follows:
    
        elevohi, /save_results, /statistics, /silent, /forMovie, /deformableFront, bgsw='HUX'
        
Up to now, three different ambient solar wind models are supported: HUX, HUXt, and EUHFORIA
For a detailed description of the ambient solar wind model maps see: 'README_AmbientSolarWinds.md'
        
