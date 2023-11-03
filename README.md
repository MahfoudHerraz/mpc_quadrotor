# mpc_quadrotor
This code implements a flatness based Model Predicitive Control algorithm for tarjectory tracking with quadrortor drones.
It is inspired from the paper: Limaverde Filho, José Oniram de A., et al. "Trajectory tracking for a quadrotor system: A flatness-based nonlinear predictive control approach." 2016 IEEE Conference on Control Applications (CCA). IEEE, 2016.

The main file (run this file) 'quadrotor_mpc.m' contains the simulation.
First, drone and simulation parameters are defined, the LTI system to be solved by MPC is defined, as well as the state obsever (Kalman filter) and mapping matrices (used in the MPC)
Second, the trajectory to be tracked is defined and variables are initialized
Third, the simulation loop is performed, solving the MPC optimization problem at each step using the fmincon Matlab function, the first control input is applied an variables are updated
Then, the total thrust, control torques and remaining states are calculated from flat outputs and their derivatives
A mix-rpm matrix allows to compute rotors speeds from total thrust and control torques

The file 'cost_function.m' implements the cost function to be minimized in the mpc optimization at each step, as well as its gradient
The file 'HessianJ.m' implements the hessian matrix of the cost function
