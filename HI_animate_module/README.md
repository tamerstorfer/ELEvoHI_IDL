# HI_animate_module

Python package for animation of Heliospheric imager modeling/prediction results

T. Amerstorfer, J. Hinterreiter, C. Möstl 

Current status (April 2020): work in progress



Install python 3.7.6 with miniconda:

Create a conda environment:

	  conda env create -f environment.yml

	  conda activate helio

	  pip install -r requirements.txt

	  


Create movie example (ffmpeg is required): 
	start python in command line

	# import the module
	from HI_animate_module import ELEvoHIEnsembleMovie as em
	# define events list (events for which the movie should be created)
    eventslist = ['20100523']
    # run the program (please make sure to set the correct paths)
    em.main(eventslist, spaceCraft='AB', scriptPath='/ELEvoHI/PredictedEvents/',
         catPath='/ELEvoHI/HI_animate_module/cats/', readData=1,
         plotBGSW=True)

	# complete list of keywords:
	#            eventsList: List with the events for which the movies should be generated
	#            spaceCraft: None or 'AB' for A and B, 'A' for A and 'B' for B
	#            readData: set to 1 if you want to create the pickle file, 0 to read already existing pickle file
	#            coordSys: HEEQ or HEE, None for HEE
	#            catPath: Path to the catalogs
	#            scriptPath: Path to the ELEvoHI ensemble output
	#            outPath: Save path for the movies
	#            plotBGSW: set to 'True' to plot the background solar wind, if available
	#            showMag: set to 'True' to plot the magnetic field legend
	#            ffmpegPath: Path to ffmpeg


Create histograms and plots for the individual events:
	start python in command line

	# import module
	from HI_animate_module import ELEvoHIVisualize_allRuns_module as ev
	# run the program with the path to the ELEvoHI runs
	ev.main('/ELEvoHI/PredictedEvents/')


