This file contains all the necessary routines to load and process data from GaitWatch and Qualisys systems.

* Content of the files:

- 'ANMS': Adaptive Nelder-Mead Simplex Algorithm. It is used to minimize the error function defined in the
          optimizer of the Kalman filter. (Auxiliary routine).

- 'calibrate_acc1D': Calibration routine for non-triaxial accelerometers. Prior to the execution of the file
	             there must exist a file in the hard drive containing the data gahtered from the calibration
                     maneuvers. (Main routine, it can be executed by the user)

- 'calibrate_acc3D': Calibration routine for triaxial accelerometers. Prior to the execution of the file there
                     must exist a file in the hard drive containing the data gahtered from the calibration
                     maneuvers. (Main routine, it can be executed by the user)

- 'calibrate_gyro1D': Calibration routine for the gyroscopes. The calibration is done uniaxially. (Main routine, 
                      it can be executed by the user)

- 'calibrate_mag3D': Calibration routine for triaxial magnetometers. Prior to the execution of the file there
                     must exist a file in the hard drive containing the data gahtered from the calibration
                     maneuvers. (Main routine, it can be executed by the user)

- 'eofKalman': Script which contains the error function which is minimized during the optimization process of 
               the Kalman filter parameters. This function is minimized using the ANMS algorithm. (Auxiliary 
               routine).

- 'gaitAnalysis': File in which the gait analysis algorithm should be written. (Main routine, it can be executed
                  by the user).

- 'GW_comm': Routine used to communicate with the GaitWatch. It can be used to change some configuration parameters
             of the sensors as well as to donwload the measured data into the computer. (Main routine, it can be 
             executed by the user).

- 'gwLibrary': File which contains all the auxiliary routines used by the main routines. 


- 'main': This is the main routine which loads the GaitWatch data, calibrates them and computes the orientation
          angles. It requires user interaction to select the data file and to select the body segments. (Main routine,
          it can be executed by the user)

- 'optimizeKF': File which calls the 'ANMS' algorithm to minimize the function in 'eofKalman'. (Auxiliary routine).

- 'param_optimization': Routine which computes the reference angles using Qualisys data. This reference is used to
                        optimize the parameters of the Kalman filter. (Main routine, it can be executed by the user).

- 'readQualisys': Routine which loads and extracts daga gathered with the Qualisys system. (Auxiliary routine).

* The 'data' folder contains:

- 'calibration' folder which contains all the calibration parameters of the sensors.
- 'raw data' folder which contains all the raw data files donwloaded from the GaitWatch, as well as the data from the
   Qualisys system measurements done the 10th and the 11th of september.
- 'workspaces' folder which contains the complete workspaces after the execution of the 'main.m' file.
- 'gwDataStruct.mat': This file contains a Matlab struct which defines the structure of the GaitWatch data.

* The 'figures' folder contains:

- 'calibratin folder' which contains a set of figures generated when executing the calibration routines.
- 'Parameter optimization' which contains a set of figures showing the pitch angle of all the leg segements for the
   experiments carried out the 11th of september.

* The 'documentation' folder contains:

- 'calibration folder' which contains a set of figures and html files which explain the calibration routines. To see
  the explanations, open the .html files.

- 'manual.pdf' A pdf file containing all the information about the GaitWatch project.


