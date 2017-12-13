Code by Rajdeep Deb:

LagrangianMainSolver.m : Final IBVP Solver. We need To specify here the mesh domain, final time when solution is desired,
                         Diffusion Coefficient and Number of time steps. Also the handle functions corresponding to 
                         Initial Condition, Boundary Condition, Source Function, Convection Vector.


SemiLagU_LFE.m : It computes the mass matrix, laplacian stiffness matrix and load vector for Iu(x-tau*v) in each of the element.


SemiLagUConstantVel_LFE.m : It computes the load vector for Iu(x-tau*v) in each of the element for stationary velocity field.


computeULagrangian.m : It computes the value of integral Iu(x-tau*v) for a given triangular element. 


computeULagrangian_ConstantVelocity.m : It computes the value of Iu(x-tau*v)*v(x) for a given triangular element. 


SourceLoad.m : Computes the source vector for a given source function


triangleSearch.m : Searches whether a given point exist inside a particular triangular element. If it exists it returns the
                  value of linear interpolated shape functions at that point.


