# PI_controller_specialcourse
This is the repository for the MATLAB code for the special course investigating PI controller instability and tuning for a floating sparbuoy wind turbine modelled as a 7 DOF system and linear rigid body dynamics.

Overall Workflow:
The main function for the code is titled "main_PI6DOF.m" and has compartmented sections of code. The general flow is inputs, then calculations, then data processing and outputs.

Controller Workflow:
The controller section includes the calculations of the required partial derivatives and also the proportional and integral constants. These values are used to run the controller in the linear system.
The linear system is run through ode4 and dqdtsparbuoy. For wind forcing only, case 3 is used to call F_wind_Region3.
F_wind_Region3 calculates the wind forcing, here the specific equation used for getting the forces can be reviewed.