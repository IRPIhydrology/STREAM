=======
STREAM- version 1.3
=======

Implementation of STREAM (SaTellite based Runoff Evaluation And Mapping) model, version 1.3 [Camici_etal_2021]_.
https://doi.org/10.5281/zenodo.4588304

In Test_data and Results you can find the dataset for testing the code and the corresponding results.

Citation
========
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4588304.svg
   :target: https://doi.org/10.5281/zenodo.4588304
   

References
==========
.. [Camici_etal_2021] Camici, S., Giuliani, G., Brocca, L., Massari, C., Tarpanelli, A., Farahani, H. H., Sneeuw, N., Restano, M., and Benveniste, J.: Synergy between satellite observations of soil moisture and water storage anomalies for global runoff estimation, Geosci. Model Dev. Discuss. [preprint], in review, 2021. https://doi.org/10.5194/gmd-2020-399


TUTORIAL for the use of the STREAM CODE
==========
STREAM is a semidistributed conceptual hydrological model for estimating runoff and river discharge from rainfall, soil moisture, air temperature and terrestrial water storage data.
The STREAM model code is distributed through M language files. It is worth to specify that the model could be run with differents interpreter of M language, like the GNU Octave (freely downloadable here https://www.gnu.org/software/octave/download). To run the code with GNU Octave it is needed to load the packages "statistics" and "optim". 
The STREAM model can be freely applied and used, please cite [Camici_etal_2021]_.
The authors are highly interested to collaborate for the understanding of the model functioning and to improve its performance and applicability.

For any questions please do not hesitate to contact:
stefania.camici@irpi.cnr.it

-----------------------------------------------------------------------------------------	

The following files are distributed:

1. M codes:

	1.1 "STREAM_semidistributed.m": STREAM model
	1.2 " cal_STREAM_semidistributed.m": code for STREAM model calibration (requires optimization toolbox)
	1.3 " run_STREAM_semidistributed.m": script for running the calibration of STREAM and creating the output figure
	1.4 " perf.m": script for computing performance scores
	1.5 " klinggupta.m": script for computing the kling-Gupta efficiency index


2. Auxiliary file:

	2.1 "Mississippi_basin.png": Mississippi sub-basins and gauging stations;
	2.2 "topology_Mississippi.csv": topology of Mississippi basin;
	2.3 "distance_Mississippi.txt": distance of each sub-basin to closure section;
	2.4 "GIUH": Geomorphological Instantaneous Unit Hydrograph;
	2.5 "X_opt_Mississippi.txt": contains the eight STREAM model parameters

	       alpha = PAR(1,i); % exponent of the infiltration 
    	        T    = PAR(2,i);  % characteristic time length
	       gamma = PAR(3,i); % parameter of GIUH
	       C     = PAR(4,i); % Celerity
               Diff  = PAR(5,i); % Diffusivity
               beta  = PAR(6,i); % coefficient relationship between slow runoff component and TWSA
               m     = PAR(7,i); % exponent relationship between slow runoff component and TWSA
               Cm    = PAR(8,i); % degree-day coefficient for snow module

3. INPUT file (example):
-----------------
	3.1 "input.mat": example file for gridded daily data for Mississippi. It contains, basin_data and temperature cells for each sub-basin.
	"basin_data" contains:
	     a) date (in numeric Matlab format)
	     b) precipitation amount (in mm)
	     c) soil moisture (dimensionless: degree of saturation)
	     d) river discharge (where availale in m3/s)
	     e) terrestrial water storage anomaly 
	"temperature" contains:
	     a) air temperature data (°C)

	3.2 "staz_checkMississippi.mat": river discharge data over the Mississippi river basin
   
STEP by STEP PROCEDURE TO RUN STREAM MODEL
==========
1. Subbasin delineation. 
	Before to run the STREAM model, a basin delineation has to be carried out. Once selected the closure sections (Nsect) over the river, this operation allows to identify 	the subbasins (Nbas) of the river basin. For instance, the basin delineation could be carried out through Qgis software (https://docs.qgis.org/3.16/en/docs/training_manual/processing/hydro.html?highlight=hydrological%20analysis).

	As output of this step, at each section selected for the analysis as well as at each subbasin should be associated a unique identification number. See for example the "Mississippi_basin.png" figure.
				 
       
2. Preparation of the input data needed for run the STREAM model. The following input data have to be created to run the STREAM model:

	2.1 "input.mat": .mat struct file of gridded daily data. It contains basin_data and temperature cells for each sub-basin.
	"basin_data" contains:
	     a) date (in numeric Matlab format)
	     b) precipitation amount (in mm)
	     c) soil moisture (dimensionless: degree of saturation)
	     d) river discharge (where available in m3/s)
	     e) terrestrial water storage anomaly 
	"temperature" contains:
	     a) air temperature data (°C)
	
       For example see the "input.mat" file provided for the Mississippi river basin.  
	

	2.2 "staz_check.mat": .mat file containing information about the river discharge data over the river basin. It contains four vectors:
	     a) Dtot:           [Nobs x 1] vector containing the date (in numeric Matlab format);
	     b) ID_bas_app:     [Nsect x 1] vector indicating the number of the basin to which each section belongs.
	     c) Q_sezcontrollo: [Nobs x Nsect] matrix containing for each section daily river discharge observations.
	     d) sez_controllo:  [Nsect x 2] matrix containing the coordinate (lon, lat) of the each section.
	
        See, for example the "staz_checkMississippi.mat" file provided for the Mississippi river basin.

	2.3 distance.txt: [Nbas x Nsect] matrix containing the distance (in km) of each subbasin to the closure sections identified over the river basin.  
			   Basins that not contribute to the closure section have a distance set equal to "-1". 

        For example see the "distance_Mississippi.txt" file provided for the Mississippi river basin.

	2.4 topology.csv: [Nbas x 6] matrix containing:
	     a) first column:  the basin sorted from the first to the last according to the numeration identified at step "Subbasin delineation";
	     b) second column: equal to the first;
	     c) third column:  connection between the basins. In other word, it specifies the number of the basin in which the river of the considered basin 
		               will continue the path.
	     d) fourth column: contains logical values to indicate if the basin is a directly draining basin (-1) or an head catchment (0).
	     d) fifth column: contains the length (km) of the river stretch belonging to the basin. 
	     d) sixth column:  contains the area (km2) of the basin.    

        See for example the "topology_Mississippi.csv" file provided for the Mississippi river basin.

3. Calibration of the model. To calibrate the model, please follow the instruction below:

	   [X_OPT]=cal_STREAM_semidistributed(input,BAS_PAR,EBRR_BASPAR,sez_outlet,bas_check,ID_bas_app)    
  
	% INPUT
	% input:   .mat struct file with input data (see above for the structure)
	% BAS_PAR:  a [3 x 1] vector containing:
		   in the 1st row the number of the subbasins (Nbas);
	           in the 2nd row the number of the section (Nsez);
	           in the 3rd row the number of the upstream input;

	% EBRR_BASPAR: [Nbas x 14] matrix containing:
		       in the 1st column the first column of topology file ;
	               from the 2nd to the Nsect+1 columns the distance as in the distance.txt;
	               in the Nsect+2 column the six column of topology file ;
	               in the Nsect+3 column the fourth column of topology file ;

	% sez_outlet: the outlet section for which to carried out the calibration model;
	% bas_check : the basin at which sez_outlet belongs;
	% ID_bas_app: vector indicating the number of the basin to which each section belongs (see above for the structure).

	% OUTPUT 
	% X_OPT: [8 x Nbas] matrix containing, for each subbasin, the calibrated model parameters.

4. Run of the model. To run the model, please follow the instruction below:
 	   
           [NS,KGE_sez,KGE_out,Qsim_out,QB_out,rr_tot]=STREAM_semidistributed(input,BAS_PAR,EBRR_BASPAR,X_OPT,sez_outlet,bas_check,ID_bas_app,FIG);


	% INPUT
	% input:   see above for the structure
	% BAS_PAR: see above for the structure
	% EBRR_BASPAR: see above for the structure
	% X_OPT: see above for the structure
	% sez_outlet: the outlet section for which to carried out the calibration model;
	% bas_check : the basin at which sez_outlet belongs;
	% ID_bas_app: vector indicating the number of the basin to which each section belongs (see above for the structure).
	% FIG: 1 for making the figure, otherwise no figure
	
	% OUTPUT

	% NS: Nash Sutcliffe Efficiency 
	% KGE_sez:  Kling Gupta Efficiency for all the Nsect sections over the basin
	% KGE_out: Kling Gupta Efficiency for the "sez_outlet" section
	% Qsim_out: Simulated total river discharge
	% QB_out: Simulated slow-flow river discharge component
	% rr_tot: Simulated gridded runoff 

	An example to load the input data, to calibrate and to run the model for the Mississippi river basin can be found within the script: "run_STREAM_semidistributed.m"   



