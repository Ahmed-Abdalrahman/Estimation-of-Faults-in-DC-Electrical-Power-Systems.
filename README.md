# PEV Charging Infrastructure Siting Based on Spatial-Temporal Traffic Flow Distribution

In this project, we propose a spatial-temporal flow capturing location model. This model determines the locations of various types of charging facilities based on the spatial-temporal distribution of traffic flows. We utilize the dynamic traffic assignment model to estimate the time-varying traffic flows on the road transportation network. Then, we cluster the traffic flow dataset into distinct categories using the Gaussian mixture model and site each type of charging facilities to capture a specific traffic pattern. We formulate our siting model as an mixed integer linear programming (MILP) optimization problem. The model is evaluated based on two benchmark transportation networks, and the simulation results demonstrate effectiveness of the proposed model.

## Files
* [Technical details](ST-FCLMpaper.pdf) -contains the technical details of the ST-FCLM model.
* [Main code](STFCLM.py) - contains the source code that interface with the traffic simulator, cluster the traffic dataset, and solve the optimization problem.
* [Results analysis](STFCLM_Figures.ipynb) - contains the analysis of the output results and data visualization.
* [Optimization model](STFCLMproblem.lp) - contains the optimization model in a way that is easier for humans to read. 
* [Traffic simulation](SUMO) - contains the traffic simulation for the Nguyen-Dupuis and the Sioux
Falls network transportation networks.

## Reference

A. Abdalrahman and W. Zhuang, “PEV charging infrastructure siting based on spatial-temporal traffic flow distribution,” IEEE Trans. Smart Grid, 2019. [Source](https://ieeexplore.ieee.org/document/8630747)

## Built With

* [SUMO](https://www.dlr.de/ts/en/desktopdefault.aspx/tabid-9883/16931_read-41000/) - Traffic flow estimation
* [Gurobi](http://www.gurobi.com/) - MILP optimization

## Authors

* **Ahmed Abdalrahman** - University of Waterloo
