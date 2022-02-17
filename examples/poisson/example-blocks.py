import flowx
import simulation


def example():
    """
    Example for block-structured Poisson Solver
    """
    nx, ny = 40, 40
    xmin, xmax = 0.0, 1.0
    ymin, ymax = 0.0, 1.0
    xblocks, yblocks = 4, 4

    # Define cell-centered variable names
    center_vars = ["asol", "ivar", "rvar", "eror"]
    poisson_vars = ["ivar", "rvar"]

    simulation_info = dict(verbose=True, poisson_solver="superlu")

    # Define boundary condition for the poisson test
    user_bc = "dirichlet"

    # Define boundary conditions for variable ivar
    bc_type = dict(ivar=[user_bc, user_bc, user_bc, user_bc])
    bc_val = dict(ivar=[0.0, 0.0, 0.0, 0.0])

    # Grid
    grid = flowx.domain.Grid(
        "cell-centered",
        center_vars,
        nx,
        ny,
        xmin,
        xmax,
        ymin,
        ymax,
        xblocks,
        yblocks,
        user_bc_type=bc_type,
        user_bc_val=bc_val,
    )

    # poisson = flowx.poisson.Poisson(grid, poisson_vars, simulation_info)

    # Compute the analytical solution
    # simulation.get_analytical(grid, 'asol', user_bc)

    # Calculate the right-hand side of the Poisson system
    # simulation.get_rhs(grid, 'rvar', user_bc)

    # Solve the Poisson system
    # ites, res = poisson.solve()

    # Compute the error (absolute value of the difference)
    # grid.compute_error('eror', 'ivar', 'asol')

    # Compute the L2-norm of the error..archive
    # l2_norm = grid.get_l2_norm('eror')
    # print(l2_norm)


if __name__ == "__main__":
    example()
