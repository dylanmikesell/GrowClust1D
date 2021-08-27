# GrowClust1D
Julia implementation of the GrowClust program for relative relocation of earthquake hypocenters based on waveform cross-correlation data. This is a testing version before later public release as a registered package. Apologies for the sparse readme; more to come.

The examples/ directory has two julia codes that run a serial version of the program without uncertainty quantification, and a multithreaded version with 100 bootstrap resamples. To run these codes, navigate to this directory and run:

] julia run_growclust-serial.jl example.serial.inp

or 

] julia -t8 run_growclust-multithread.jl example.multithread.inp

In the second example, I specified the usage of 8 threads (you can decide based on your resources).