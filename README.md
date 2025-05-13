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

1. Folder List

   A. registration: Contains MATLAB scripts for registering 2-photon data to Allen Institute CCF

   B. python: Contains Python and Jupyter Notebook scripts for pre-processing the data

   C. utils: Contains helper functions used throughout the analysis

   D. ridgeModel: Contains the git repository from https://github.com/musall/ridgeModel


---------------------------
METHODOLOGICAL INFORMATION
---------------------------

1. Description of methods used for collection/generation of data: 

	Behavior imaging

		Mice were head-restrained in a plexiglass apparatus and the scene was illuminated with 850nm LEDs (Gupta and Murphy 2025).
		A monochrome camera (Omron-Sentech, STC-MBS43U3V) connected to a NVIDIA Jetson was equipped with a bandpass infrared
		filter (840-865 nm, Bock Optronics) and positioned directly in front of the animal to record 8-bit behavior videos at 90
		frames per second with a resolution of 640 x 320 pixels (Figure 1A). 
	
	1-photon imaging

		Dorsal cortical GCaMP activity was recorded using similar methods as described previously (Ramandi et al. 2023). Briefly,
		images of the cortical surface (Figure 3A) were projected through a pair of back-to-back photographic lenses
		(50 mm, 1.4 f:135 mm, 2.8 f or 50 mm, 1.4 f:30 mm, 2 f) onto a 1M60 Pantera CCD camera (Dalsa). GCaMP was excited with a
		blue LED (Luxeon, 473 nm) with a band-pass filter (Chroma, 467–499 nm) delivered to the surface of the cortex through a
		liquid light guide (Thorlabs). GCaMP fluorescence emission was filtered using a 510–550 nm band-pass filter (Chroma).
		12-bit images were collected at 30 frames per second using XCAP imaging software using 8x8 pixel on-chip binning, yielding
   		images of size 128x128 pixels. 

	2-photon imaging

		Two-photon imaging was performed using the Diesel2p microscope (Yu et al. 2021) controlled with MATLAB software Scanimage
		(Vidrio). Two-photon laser excitation was provided by an 80-MHz Newport Spectra-Physics Mai Tai HP 1020. The Diesel2p’s
		dual scan engines were utilized to maximize the imaging area and imaging regions of interest (ROIs) were chosen across mice
		and imaging sessions to cover as many cortical regions as possible. Each scan path consisted of a resonant scanner
		(CRS 8 KHz, Cambridge Technologies), and two XY galvanometers (8320K, Cambridge Technologies). While imaging area sizes
		sometimes varied between scan engines, the spatial resolution was held constant across scan arms, resulting in acquisition
		rates that varied, with larger areas collected at lower frames per second. Photon signal was collected by photomultiplier
		tube (H11706P-40, Hamamatsu) and amplified with a high-speed current amplifier (HCA-400M-5K-C, Laser Components). Imaging
		was performed with 920 nm laser at ~80 mW excitation out of the front of the air objective (0.55 NA). 

3. Methods for processing the data: 

	Identification of grooming behaviors

		Grooming behaviors were identified using two independent approaches - manual labelling (Figure 1-1A) and unsupervised
		clustering (Figure 1D). Final results are presented using the unsupervised clustering approach, however both labelling
		strategies demonstrated substantial agreement (Figure 1D, Figure 1-1B) and the results presented throughout the paper were
		consistent across approaches. For manual labelling, the start and end index for each grooming event were saved and manually
		curated using Behavioral Observation Research Interactive Software (BORIS) (Friard and Gamba 2016). Grooming component
		behaviors were classified by visual inspection as one of seven distinct behaviors, including unilateral facial strokes:
		right and left; and bilateral facial strokes: elliptical, elliptical asymmetric, right asymmetric, left asymmetric, and
		large bilateral (Figure 1-1B). BORIS was also used to label licking behaviors where the tongue clearly protruded from the
		mouth, as well as the precise timepoints in which the full weight of the drop makes contact with the mouse.

		For unsupervised labelling, the position of the paws and nose were tracked using DeepLabCut (Mathis et al. 2018). A single
		point was used for each paw and around 20-60 frames from each trial and mouse were selected for labeling. The network was
		trained with default parameters for 1030000 iterations. The start and end frame index of each grooming behavior obtained
		using BORIS were used to extract paw position and velocity information during each behavior event. This information was
		distilled down to  44 features for each behavior event (Table 1), yielding a matrix of size Nx44, where N is the number of
		grooming behaviors exhibited across all mice. This feature matrix was then reduced to 2 dimensions using Uniform Manifold
		Approximation and Projection (UMAP) for dimension reduction. The resulting embedding exhibited several spatially distinct
		clusters, which were extracted with hierarchical clustering (Figure 1D). These clusters correspond with distinct and
		stereotyped patterns of paw movements (Figure 1E, Figure 1-1C).

	1-photon image pre-processing
   
		Single photon wide-field fluorescence data was compressed using singular value decomposition (SVD), which yielded U, the
		matrix of pixels × components; V, the matrix of components x time; and s the diagonal matrix of singular values. The top
		1,000 components were retained, and s was multiplied into V, which was then upsampled to 90 fps to match the behavior
   		video. All subsequent analyses, such as temporal filtering and ridge regression, were performed on tThe matrix s*V was
   		then high-pass filtered at 0.01 Hz using a second-order zero-phase Butterworth filter. Wide-field fluorescence data were
   		combined within mice across separate recordings, by recasting experiment-specific SVD components into a master SVD basis
   		set (Peters et al. 2021). To create the master SVD basis set, spatial components from each imaging session were first
   		aligned to a reference session and concatenated, then SVD was performed on the resulting matrix and the top 1,000
   		components were retained. Temporal components for each experiment then were recast into the master basis set as described
   		in (Peters et al. 2021). Wide-field fluorescence data for each mouse were then aligned to the Allen Institute Common
   		Coordinate Framework (Q. Wang et al. 2020) using a control-point based registration method with key points placed along
   		the superior sagittal sinus, bregma, and the space between olfactory bulbs and frontal cortex (Saxena et al. 2020). A
   		mask was drawn manually over the cortex of the reference image for each mouse to discard pixels outside of the cortex
   		from analysis.

	2-photon image processing
   
		2-photon images were spatially binned by a factor of 2x2 to improve signal to noise ratio and reduce file sizes for
   		subsequent operations. After binning, neuron diameters were typically ~5-8 pixels. Motion correction and neuron detection
   		were then performed using Suite2p. Neuronal regions of interest were curated after visual inspection of their shape and
   		fluorescence signals. Neuronal fluorescence signals were consolidated by upsampling the signals obtained with the lowest
   		frame rate to match those acquired with the maximum frame rate. To relate neuronal activity to behavior, the consolidated
   		fluorescence signals were upsampled again to match the framerate of the behavior video.


5. Instrument- or software-specific information needed to interpret the data: 

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


6. Environmental/experimental conditions: 

		After an initial baseline period (30s for 1-photon imaging or 5 mins for 2-photon imaging), grooming behavior
		was evoked by delivering a drop of water onto the orofacial area once per minute for 20 minutes. Each water drop stimulus
		was preceded by a 1-second 10kHz tone followed by a 1-second delay period (Figure 1B). For each mouse, at least two
		spontaneous sessions, which contained the audio cue without the water drop stimulus, were recorded before the evoked
		grooming sessions (Figure 1C).

8. Describe any quality-assurance procedures performed on the data: 

	Synchronization between devices

		Video acquisition was triggered by the Jetson which also started behavior video recording, and after a few seconds of
		delay, the blue cortical excitation LED and the infrared behavior LED were simultaneously turned on by TTL from the
		Jetson. At the end of the trial, the cortical and behavior LEDs were turned off simultaneously prior to stopping
		acquisition. Frames were synchronized across behavior and brain cameras by matching the illuminated frames at the
		start and end of the trial.

   		Behavior video acquisition and resonant scanning were initiated prior to collecting 2-photon images. Acquisition of 2-photon
   		images was triggered using TTL output from the Jetson which simultaneously turned on infrared LEDs which illuminated the
   		scene, and at the end of the trial, 2-photon image acquisition was terminated by TTL and the infrared LEDs were
   		simultaneously turned off. Behavior video acquisition ended following a short delay after turning off the LEDs and 2-photon
   		images were synchronized to the behavior video by matching to the illuminated frames.

   	Alignment of 2-photon images to Allen Common Coordinate Framework

		Prior to 2-photon imaging experiments, the locations of the forelimb and hindlimb somatosensory cortex were determined by
		sensory mapping. Sensory mapping experiments were conducted under similar conditions as the single-photon imaging
		experiments, only mice were anesthetized with 1.2% isoflurane in oxygen and images were acquired at 40 fps (Figure 4-1D).
		Piezoelectric stimulators touching the forelimb or hindlimb were used to record somatosensory-evoked stimulus responses
		(Y. Xie et al. 2016). Averages of responses to sensory stimulation were calculated from 20 trials of stimulation with an
		interstimulus interval of 10 s (Figure 4-1E). After each 2-photon imaging session, a 5x5mm2 image of the cortical surface
		was taken with the Diesel2p’s linear scanners. ROIs captured with the resonant scanners were registered to this larger
		field-of-view 2-photon image, which was then registered to a 1024x1024 pixel template image (8.2x8.2mm2), obtained on the
		single-photon imaging system (Figure 4-1A-C). The single-photon template image was then registered to the dorsal projection
		of the Allen Common Coordinate Framework using control points placed in the forelimb and hindlimb somatosensory cortex
		(Figure 4-1F). Each image registration process yields a transformation matrix, which was used to transform each neuron’s
		pixel location in the resonant image obtained from Suite2p, to its location in the Allen Common Coordinate Framework.

10. People involved with sample collection, processing, analysis and/or submission: 
Nicholas Michelson

--------------------------
DATA-SPECIFIC INFORMATION 
--------------------------

Each experiment includes:


Note: 

