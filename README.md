# RNA-transcription-modelling-and-DS
Model Implementation and Design Space Construction to support 'Quality by Design modelling for rapid RNA vaccine production against emerging infectious diseases'

This README document explains the function of each item within this repository and how the items are linked to each other

### Dependencies
- *environment.yml*: Lists all packages plus the version used. Using `conda env create --file environment.yml` in the anaconda prompt will recreate the author's dependencies

- *requirements.txt*: Lists all packages plus the version used. If the reader prefers pip, a virtual environment should be created and the author's dependencies recreated with `pip install -r requirements.txt`

### Documentation
- *Documentation.md*: Gives scientific background reading on the functions and main scripts, and comments on some of the code implementation

### Function scripts
-  *mass_balances.py*: Includes the solution component concentration mass balances that need to be solved for either constant pH or overall H concentration at each time step within *odes_and_curve_fitting_functions*

- *odes_and_curve_fitting_functions.py*: Includes all functions from the definition of the system of ODES definining the kinetic expressions of the overall component species, to the functions looking for an initial guess to solve the mass balance considerations at each time step, to the solution of the ODEs for a given data sample. The functions in this script are called upon whenever simulations are required.

- *data_functions.py*: Includes the function that loads the data from the *Experimental_Data.xlsx* file and extracts relevant variables, and the function that plots the RNA yield for each experimental sample

### Main scripts
It is recommended to run the scripts in the following order to replicate the workflow followed by the authors, but each of these scripts can be run independently. The scripts that are marked as optional do not contain figures that directly support the manuscript.
1. *curve_fitting_exploration_prediction.py*: Uses the *scipy.optimize.curve_fit()* tool to estimate the model parameters. It gives the optimal mdoel parameters, their standard deviation and correlations. It then plots the simulated versus experimental results for each data sample, produces a prediction error plot, and shows the explored input space, namely the data samples in the process parameter space.

2. *optimal_parameters_exploration.py* (optional): Explores in more depth the dependency of the RNA yield on the process parameters, namely initial Magnesium, T7RNAP and NTP concentrations at the curve_fitting's optimal model parameter values.

3. *cross_validation.py* (optional): Does 10-fold cross validation. It randomly splits the experimental samples into 10 folds. It performs parameter estimation (curve_fit()) 10 times using 9 folds as the training set and leaving the remaining set to test the fit. It then plots the simulated versus experimental results for each data sample as well as the prediction error plot

4. *Cost_Yield_Plots.py*: Fixes the Mg concentration and produces a 3D figure showing the cost per gram of RNA as a function of T7RNAP and Mg concentration to find the cost-optimal operating point and its associated cost. This graph is colour-coded according to the abolute RNA yield

5. *3D_deterministic_DS.py*: Creates a 3D grid of the input process parameter space consisting of initial Mg, NTP and T7RNAP concentrations and shows the design space, the grid points that meet a certain yield threshold at the previously found optimal model parameters. This figure is colour-coded according to the absolute RNA yield.

6. *2D_probabilistic_DS.py*: Fixes T7RNAP and performs 50 Monte Carlo simulations where model parameters are sampled with 20% standard deviation around their experimental optimum to get the probability of a point on the Mg-NTP grid reaching a required yield threshold.

### Excel files
- *Experimental_Data.xlsx*: Includes all the data samples with the relevant conditions used for parameter estimation. These are part of a bigger Design of Experiments dataset by Prof. Robin Shattock's group. This file is loaded within the *data_functions.py* script, which itself is called in *cross_validation.py*, *curve_fitting_exploration_prediction.py* and *optimal_parameters_exploration.py*

- *3d_DesSpace_new.csv*: Includes the results of a previously simulated RNA yield 3d-deterministic design space. These can be loaded within *3D_deterministic_DS.py* if only the plots need to be generated. These can also be overwritten in the same file with simulation results from different parameters

- *y_sim_prob_Design.csv*: Includes the results of 50 previous Monte Carlo simulations with 20% error on the model parameters. These are not directly used in any function but might come in handy to generate the equivalent of *prob_Design.csv* should the latter become inapplicable if the code is altered.These can be loaded in *2D_probabilistic_DS.py* and overwritten in the same file with simulation results from different parameters

- *prob_Design.csv*: Includes the probabilities of the sample grid points from 50 previous Monte Carlo simulations with 20% error on the model parameters. These can be loaded within *2D_probabilistic_DS.py* if only the plot needs to be generated. These can also be overwritten in the same file with simulation results from different parameters

### To ignore
- *pycache* folder: Can be ignored or deleted, as this is generated when running the scripts the first time in Python
