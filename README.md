# Estimation of Faults in DC Electrical Power Systems.

Course project: ECE 602- Winter 2017

The goal of this project is to replicate and analyze the estimation methods provided in [1]. It is noteworthy to say that the results provided in [1] was experimentally achieved, and do to lack of equipment we solely focus on simulation studies.

The main goal of this paper is to estimate faults in a DC electric circuit via optimization methods. A DC electric circuit consists of sources, loads and switching elements. voltage and current measurements are implemented on some specific nodes in the system. 

## Files
* [Technical details](paper.pdf) -contains the technical details of the model.
* [Main code](main.py) - contains the source code that solve the optimization problem.
* [Results analysis](main.html) - contains the analysis of the output results and data visualization.

## Reference

[1] D. Gorinevsky, S. Boyd and S. Poll, "Estimation of faults in DC electrical power system," 2009 American Control Conference, St. Louis, MO, 2009, pp. 4334-4339.

## Built With

* [MATLAB]  
* [CVX](http://cvxr.com/cvx/) - Convex optimization

## Authors

* **Ahmed Abdalrahman and Farid Farmani** - University of Waterloo
