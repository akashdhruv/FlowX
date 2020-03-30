"""Routine to solve the Poisson system with iterative or direct solvers."""
  # (In that case I calculate the RHS and the analytical solution here, because there is no need of user input,
  # but in the future it will be defined by the user in the jupyter script)

import numpy
import flowx
import time
import examples.poisson.simulation as simulation 

def solver(grid, asol, ivar, rvar, user_method, user_bc, error,  maxiter=3000, tol=1e-9, verbose=False):
    """Solve the Poisson system using iterative or direct solvers.
    
    Arguments
    ---------
    grid : Grid object
        Grid containing data.
    asol : string
        Name of the grid variable of the analytical solution.
    ivar : string
        Name of the grid variable of the numerical solution.
    rvar : string
        Name of the grid variable of the right-hand side.
    user_method : string
        Name of the numerical method to solve the Poisson system.
    error : string
         The difference between the analytical and the numerical solution at every iteration. 
    maxiter : integer, optional
        Maximum number of iterations;
        default: 3000
    tol : float, optional
        Exit-criterion tolerance;
        default: 1e-9

    Returns
    -------
    iters: integer
        Number of iterations computed.
    difference: float
        Final difference between u at time step n+1 and n.
    verbose : bool, optional
        Set True to display convergence information;
        default: False.
    error_max : float
        The maximum difference between the analytical and the numerical solution at every step
    duration : float
        Time of solution when completed  
    number_iter : numpy.ndarray
        1D array of integers. 

    """

    #get analytical solution for u
    simulation.get_analytical(grid, asol, user_bc)

    if (user_method == 'jacobi'):
 
        #initialise u
        u = grid.get_values(ivar)

        #set and get RHS
        simulation.get_rhs(grid, rvar, user_bc)
        b = grid.get_values(rvar)

        #get dx, dy
        dx, dy = grid.dx, grid.dy

        #initialise variables for the loop
        iters = 0
        difference = tol + 1.0
        error_max = numpy.zeros(maxiter)
        number_iter = numpy.zeros(maxiter)
    
        start = time.time()
        #Solve the system
        while iters < maxiter and difference > tol:
            u_old = numpy.copy(u)  #previous solution
            #calcualte u numerical
            u[1:-1, 1:-1] = (((u_old[1:-1, :-2] +
                                 u_old[1:-1, 2:]) * dy**2 +
                                (u_old[:-2, 1:-1] +
                                u_old[2:, 1:-1]) * dx**2 -
                                b[1:-1, 1:-1] * dx**2 * dy**2) /
                               (2 * (dx**2 + dy**2)))

            #add BC
            grid.fill_guard_cells(ivar)

            #difference
            difference = (numpy.sqrt(numpy.sum((u - u_old)**2) /
                    ((grid.nx + 2) * (grid.ny + 2)))) 
           
            #error
            grid.get_error(error, ivar, asol)
            error_max[iters] = grid.get_linfinity_norm(error)
            iters += 1
            number_iter[iters - 1] = iters

        end = time.time()
        duration = end-start
        
        #print
        if verbose:
            print('Jacobi method:')
            if iters == maxiter:
                print('Warning: maximum number of iterations reached!')
            print('- Number of iterations: {}'.format(iters))
            print('- Final difference: {}'.format(difference))
            print('- Time: {}'.format(duration))

        return iters, duration, error_max, number_iter

    elif (user_method == 'gauss_seidel'):

        #initialise u
        u = grid.get_values(ivar)

        #set and get RHS
        simulation.get_rhs(grid, rvar, user_bc)
        b = grid.get_values(rvar)
        
        #get dx, dy
        dx, dy = grid.dx, grid.dy

        #initialise variables for the loop
        iters = 0
        difference = tol + 1
        error_max = numpy.zeros(maxiter)
        number_iter = numpy.zeros(maxiter)

        start=time.time()
        #j:column, i:row
        while iters < maxiter and difference > tol:
            u_old = numpy.copy(u)  # previous solution
            for j in range(1, grid.ny + 1):
                for i in range(1, grid.nx + 1):
                    u[i, j] = (((u[i, j - 1] + u[i, j + 1]) * dy**2 +
                                  (u[i - 1, j] + u[i + 1, j]) * dx**2 -
                                  b[i, j] * dx**2 * dy**2) /
                                 (2 * (dx**2 + dy**2)))

            grid.fill_guard_cells(ivar)

            difference = (numpy.sqrt(numpy.sum((u - u_old)**2) / ((grid.nx + 2) * (grid.ny + 2))))

            #error
            grid.get_error(error, ivar, asol)
            error_max[iters] = grid.get_linfinity_norm(error)
            iters += 1
            number_iter[iters - 1] = iters        

        end=time.time()
        duration=end-start
        
        if verbose:
            print('Gauss-Seidel method:')
            if iters == maxiter:
                print('Warning: maximum number of iterations reached!')
            print('- Number of iterations: {}'.format(iters))
            print('- Final residual: {}'.format(difference))
            print('-Time: {}'.format(duration))

        return iters, duration, error_max, number_iter

    elif (user_method == 'direct_inversion'):
 
        from scipy.sparse.linalg import spsolve

        #set and get RHS
        simulation.get_rhs(grid, rvar, user_bc)
        b = grid.get_values(rvar)*grid.dx**2
        b = b[1:-1,1:-1].flatten()

        #get A matrix
        start_m = time.time()
        A, D = flowx.poisson.construct_matrix(grid, user_bc)
        end_m = time.time()
        duration_m = end_m - start_m
 
        #solve the system
        start = time.time()
        u = spsolve(A, b)
        end = time.time()
        duration = end - start
        
        #map u to the 2D grid
        u = numpy.reshape(u,(grid.nx,grid.ny))
        temp = grid.get_values(ivar)
        temp[1:-1,1:-1] = u
        grid.fill_guard_cells(ivar)    
 
        #print
        if verbose:
            print('Direct Inversion Method:')
            print('- Time of solver: {}'.format(duration))
            print('- Time of matrix construction: {}'.format(duration_m))   

        return duration

    elif (user_method == 'super_lu'):
        
        from scipy.sparse.linalg import splu

        #set and get RHS
        simulation.get_rhs(grid, rvar, user_bc)
        b = grid.get_values(rvar)*grid.dx**2
        b = b[1:-1,1:-1].flatten()

        #get A matrix
        start_m = time.time()
        A, D = flowx.poisson.construct_matrix(grid, user_bc)
        end_m = time.time()
        duration_m = end_m - start_m

        #solve the system
        start = time.time()
        LU = splu(A)
        u = LU.solve(b)
        end = time.time()
        duration = end - start
        
        #map u to the 2D grid
        u = numpy.reshape(u,(grid.nx,grid.ny))
        temp = grid.get_values(ivar)
        temp[1:-1,1:-1] = u
        grid.fill_guard_cells(ivar)

        #print
        if verbose:
            print('Solution using SuperLU factorization:')
            print('- Time: {}'.format(duration)) 
            print('- Time of matrix construction: {}'.format(duration_m))
   
        return duration

    elif (user_method == 'conjugate_gradient'):

        from scipy.sparse.linalg import cg

        #set and get RHS
        simulation.get_rhs(grid, rvar, user_bc)
        b = grid.get_values(rvar)*grid.dx**2
        b = b[1:-1,1:-1].flatten()

        #get A matrix
        start_m = time.time()
        A, D = flowx.poisson.construct_matrix(grid, user_bc)
        end_m = time.time()
        duration_m = end_m - start_m

        #solve the system
        iters = 0
        def callback(xk):
            nonlocal iters
            iters += 1
        start = time.time()
        u = cg(A, b, maxiter=maxiter, tol=tol, callback=callback) 
        end = time.time()
        duration = end - start
        u = u[0:-1]
        
        #map u to the 2D grid
        u = numpy.reshape(u,(grid.nx,grid.ny))
        temp = grid.get_values(ivar)
        temp[1:-1,1:-1] = u
        grid.fill_guard_cells(ivar)

        #print
        if verbose:
            print('Solution using Conjugate Gradient method:')
            print('- Time: {}'.format(duration))
            print('- Time of matrix construction: {}'.format(duration_m))

        return duration
