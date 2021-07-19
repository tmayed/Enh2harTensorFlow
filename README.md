# Enh2harTensorFlow

Machine learning optimisation using TensorFlow for enhanced second harmonic grating parameters. Parameter data is generated in enh_rand_sim.py and outputed in a csv. Simulation solvers are written in Fortran (enh_solvers.f90) and compiled using f2py. Use Makefile to compile. TensorFlow model is trained in enh_tf.ipynb from csv data file.
