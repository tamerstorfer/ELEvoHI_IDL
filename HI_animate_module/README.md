# HI_animate_module

Python package for animation of Heliospheric imager modeling/prediction results

by C. Möstl, T. Amerstorfer, contributors: J. Hinterreiter

Current status (April 2020): work in progress



Install python 3.7.6 with miniconda:

Create a conda environment:

	  conda env create -f environment.yml

	  conda activate helio

	  pip install -r requirements.txt

	  



To create the movie change in 'ELEvoHIEnsembleMovie.py' the parameters in 'main'.
Then run: 'python ELEvoHIEnsembleMovie.py'


To create the Histograms and other images change the path to the runs in 'ELEvoHIVisualize_allRuns_module.py'
Then run: 'python ELEvoHIVisualize_allRuns_module.py'

