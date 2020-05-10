Program Structure
==================

* main.c : sets the testcases. Once all the structure are created, the time integration scheme can begin. Its structure is the following :

      1. Update cells
      2. Update neighbors
      3. Compute pressure according a state equation
      4. Time integration
      5. Update the animation
      6.  If t < tEnd : go to 1.

* time_integration.c : provides all the useful functions to compute the time integration and the calculation of the pressure field.
                   The time integration follows the scheme :

                   1. Update density
                      1a. Density correction based on CSPM
                   2. Update velocity
                      2a. if requires the pressure gradient is corrected using CSPM.
                      2b. XSPH is also available by the set of parameter eta, set to 0 for zero effect and to 0.5 for commun value.
                   3. Update position
                   4. Boundary conditions
                   5. show animation

* particle.c : this code describes the structure of a particle in this framework. It also provides all the necessary functions to compute
           the particle neighbords. The particle structure describes the physics information that each particle weights.

* boundary.c : imposes Dirichlet condition by imposing a boundary velocity. For now, only the tangent velocity is updated but an update of the normal velocity is straightforward.

* SPH_operator.c : computes the gradient, laplacian and divergence based on SPH method. The kernel is defined in kernel.c.

* kernel.c : returns the value of the kernel function and its derivative used in the SPH approximation.

* animation.c : create the window to see your results in real time. It is possible to save your results by decommenting lines 389-391 and 414.

* vector and list are codes to abstract some expression.
