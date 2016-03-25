About generator: Generates simulation to text file.

TO ALLOW THIRD DIMENSION OF PARTICLES: Uncomment lines in function generateRandomTestParticles() in generator.cpp.

How to compile: make compile


How to run the generator:

1) ./generator

Result: File "output.txt" with default settings of simulation.

2) ./generator output_file.txt

Result: The generator will interactively ask for simulation parameters and it will save the resulting simulation in "output_file.txt".


3) ./generator amount sim_steps width height depth max_speed max_weight output_file.txt

Result: The generator will create a file "output_file.txt" with simulation settings according to command line arguments.



DEFAULT SIMULATION PARAMETERS (specified in initDefaultSimConfig() in generator.cpp (maybe move it entirely to constructor of SimConfig?)):

Amount of particles = 100

Simulation steps = 1000

Width = 640

Height = 640

Depth = 500

Max speed = 4

Max weight = 65536