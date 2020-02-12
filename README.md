# LinearMPC
implementation of model predictive control for a linear system.

The plant model used in the simulation code is the quadruple-tank process as a MIMO linear system [1].

For a quadratic programming solver, this sample code employes 'quadprog' included in MATLAB optimization toolbox, but any other optimization solvers for linear or convex problem can be configured for this program. 

## Simulation Result
![alt text](https://github.com/jinsung20/LinearMPC/blob/master/Result.png)

## Requirements
MATLAB Optimization Toolbox
or
3rd party quadtratic programming solver, etc.

## Reference
[1] K. H. Johansson, "The quadruple-tank process: a multivariable laboratory process with an adjustable zero," IEEE Transactions on Control Systems Technology, vol. 8, no. 3, pp. 456-465, May 2000.
