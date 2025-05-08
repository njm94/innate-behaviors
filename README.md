--------------------
GENERAL INFORMATION
--------------------

1. Title of Dataset:  

2. Author Information
	A. Principal Investigator Contact Information
		Name: Timothy Murphy
		Institution: University of British Colombia
		Email: thmurphy@mail.ubc.ca

	B. Associate or Co-investigator Contact Information
		Name: Nicholas Michelson
		Institution: University of British Colombia
                Email: njm89@student.ubc.ca


3. Date of data collection: 2023-11-06 through 2024-08-02

4. Information about funding sources that supported the collection of the data: Canadian Institutes of Health Research (CIHR), Brain Canada, Canadian Open Neuroscience Platform 


---------------------------
SHARING/ACCESS INFORMATION
---------------------------

1. Licenses/restrictions placed on the data: 
These data are available under a Creative Commons Attribution 4.0 International (CC BY 4.0) license <https://creativecommons.org/licenses/by/4.0/> 

2. Links/relationships to ancillary data sets or software packages: 



3. Recommended citation for this dataset: 



---------------------
DATA & FILE OVERVIEW
---------------------

1. File List

   A. Folder:     
      Short description:       

   B. Folder:     
      Short description: 

   C. Folder name:      
      Short description: 
	Note: 

   D. Folder name: 
      Short description: 
	Note: 


---------------------------
METHODOLOGICAL INFORMATION
---------------------------

1. Description of methods used for collection/generation of data: 

	Behavior imaging

		Mice were head-restrained in a plexiglass apparatus and the scene was illuminated with 850nm LEDs (Gupta and Murphy 2025).
		A monochrome camera (Omron-Sentech, STC-MBS43U3V) connected to a NVIDIA Jetson was equipped with a bandpass infrared
		filter (840-865 nm, Bock Optronics) and positioned directly in front of the animal to record 8-bit behavior videos at 90
		frames per second with a resolution of 640 x 320 pixels (Figure 1A). After an initial baseline period, grooming behavior
		was evoked by delivering a drop of water onto the orofacial area once per minute for 20 minutes. Each water drop stimulus
		was preceded by a 1-second 10kHz tone followed by a 1-second delay period (Figure 1B). For each mouse, at least two
		spontaneous sessions, which contained the audio cue without the water drop stimulus, were recorded before the evoked
		grooming sessions (Figure 1C). 

3. Methods for processing the data: 


4. Instrument- or software-specific information needed to interpret the data: 

	Instruments

	Software

		A. Data preprocessing
			- Python		3.7 
			- numpy 		1.16.4
			- opencv-python		4.1.0.25
			- scipy 		1.2.1
			- roipoly 		0.5.0
			- scikit-learn 		0.21.2
			- matplotlib 		3.1.0

		B. Data analysis


5. Environmental/experimental conditions: 


6. Describe any quality-assurance procedures performed on the data: 


7. People involved with sample collection, processing, analysis and/or submission: 
Nicholas Michelson

--------------------------
DATA-SPECIFIC INFORMATION 
--------------------------

Each experiment includes:


Note: 

