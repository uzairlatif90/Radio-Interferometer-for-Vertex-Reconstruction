# Interferometer for Radio Detectors

In order to do vertex position reconstruction we need to carry out interferometry. Interferometery is basically a minimization problem where have a pulse arrival times from data which we need to match with pulse arrival times for a given transmitter position using simulation (i.e., raytracing).

The minimization happens over the spherical coordinates (theta,phi,R) of the transmitter position. The coordinates are scanned over and raytracing is using to calculate the pulse hit times for receiving antennas for each transmitter position.

The process ends when the minimizer finds the transmitter position that gives the closest match to the pulse hit times. The Tx coordinates that give the best match are then taken to be the  reconstructed coordinates of the original vertex.

## A quick overview

The minimization over Theta, Phi and R is broken down into:

- A 1-D minimization for R and then
- A 2-D minimization over (θ,Φ) for each value of R
- An initial guess value of theta, phi and R is passed onto the minimizer to start the process. The radial limit is set to be +-500 m around the guess distance value.
- Currently I am just assuming that we perfectly get D and R times every time a ray hits an antenna.

I have set my origin to be at the station center. The center is basically at the average coordinate of all the receiving antennas. So theta<90 deg is above the station and theta>90 deg is below the station. Theta 0 deg is straight up and theta 180 deg is straight down. Distance is always given in meters.

In all of this I use my own analytic raytracing which can be accessed [here](https://https://github.com/uzairlatif90/IceRayTracing). By default I use the AraSim South Pole exponential refractive index profile which is given by n(z)=A+Be^{Cz} where A=1.78, B=-0.43, C=-0.0132 m^-1.

### Minimization function

I minimize (or carry out the interferometry) the hit times of all the recieving channels/antennas by calculating the following pseudo-chisquared value:

<img src="https://render.githubusercontent.com/render/math?math=\chi^2=\frac{\sum^{N}_{i} \left(\frac{\sum^{N}_{j} (dt(i,j)_{data}-dt(i,j)_{sim})^2 w_{j}}{\sum^{N}_{j}w_j} \right)w_i}{\sum^{N}_{i}w_i}\text{ here }dt(i,j)=T_i-T_j,~~w_i=SNR_i">

If we also include R (or the second pulse) then the variables will be updated to:

<img src="https://render.githubusercontent.com/render/math?math=N &=N_D+N_R\\ \sum_i^N w_i &= \sum_i^{N_{D}} w_{D,i} + \sum_i^{N_R} w_{R,i}">

## Scripts

- RunInterferometer.C : This script takes in as arguments of a given vertex and then tries to reconstruct the given position using interferometry. In this script the user has the option of scrambling hit times and changing the ice model (i.e., have a different a ice model for generating hit times and a different one for reconstructing the position). By default an ARA type station is chosen for this exercise though the user can define their own desired radion station in AntennaConfig.hh and then work with that. AntennaConfig.hh contains pre-defined station layout for ARA and RNO-G stations.  

- RunInterferometerMultTx.C : This script does the same thing as RunInterferometer.C but for multiple vertex position. The user can define a regular grid in 3-D where multiple vertex positions can be simulated and then the interferometer will try to reconstruct them. The grid spacings and dimension can be altered by the user. By default an ARA type station is chosen for this exercise though the user can define their own desired radion station in AntennaConfig.hh and then work with that. AntennaConfig.hh contains pre-defined station layout for ARA and RNO-G stations.

- RunInterferometerARA.C : This script tries to reconstruct an ARA CP with an ARA station. It uses antenna coordinates obtained from AraRoot. All the hit times in the recieving channels are generated using the raytracer. In this script the user has the option of scrambling hit times and changing the ice model (i.e., have a different a ice model for generating hit times and a different one for reconstructing the position).

- ClusterScripts : This folder contains scripts which are just implementations of the interferometer for some specific cases like for looking ARA data and ARA SPICE data. If the user wants to use to the interferometer to simulate several vertex positions than the SIM folder provides some scripts which can be run on a computing cluster to reconstruct multiple vertex positions at the same time.

## Prerequisites
You will need to have a functioning installation of [GSL](https://www.gnu.org/software/gsl/) ([2.6](https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz) is verified to work).
- You will need to set the environment variable `GSLDIR` to your local installation of GSL.
- You will also need to have `GSLDIR` in your `LD_LIBRARY_PATH`.
- For Mac users: you can locate your GSL installation via `gsl-config --prefix`, i.e. `export GSLDIR=$(gsl-config --prefix)`
- If you have Ubuntu or Linux you can skip all of the above and get it from the repository by doing: `sudo apt install libgsl-dev`

## How to install and run the scripts

### RunInterferometer.C as a standalone package
To run, you have to do:
- `root -l RunInterferometer.C`

### RunInterferometerMultTx.C as a standalone package
To run, you have to do:
- `root -l RunInterferometerMultTx.C`

### RunInterferometerARA.C as a standalone package
To run, you have to do:
- `root -l RunInterferometerARA.C`
