# OWC-Internship
In this repertory, You will find 5 folders in which there are programs that I used to have results during my Internship in Tecnalia. See the report for more detail. 
You will find in the following folders : 

- Air compressibility : It gives you a spring-like effect of air compressibility introduced by a phase shift between
pressure and flow rate. It plots also all the models (incompressible, isentropic, and polytropic) on the same figure. Relative errors between models are calculated. You can do it for the Biradial and the Wells. For a type of turbine we have in input a certain type of wave elevation depending on the sea state that you wish. The Motion of the water surface column comes from a Simulink block which can be found in thermo_validity and Validity_hydrodynamic and it is a 100% simulation. 

- KnowTheSwell : It gives the heigth, period (and even more information but it is not used) of the waves in the sea every hours for the month of February, march and April. 
Err_Allfev.mat are results from Thermo_validity, it is interesting to look at the swell that we had to understand the results.

- Thermo_validity : We skip the modelisation of the Hydrodynamic to evaluate the validity of the thermodynamic models that were created, by giving in input post processed measurements that were collected in Open sea. We also compare the model with other collected data from Mutriku's chamber to have the relative pressure in it. Relative errors between the data and the models (isen and poly) are computed to create the trend (see report) for the month of February. 

- FetchDataBase : A huge data base is accessed via Internet(free access but controlled). Smaller matrixes are made to work faster.

- Hydrodynamic Model : The Hydrodynamic model is compared with collected data. We give in input to the model measurements from the Swell in open sea and we compare his output (motion which incudes position of the water in the OWC and its velocity) with the motions that we record in the OWC. 
