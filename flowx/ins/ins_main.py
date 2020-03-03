"""Module for incompressible Navier Stokes equations"""

from flowx.ins.ins_interface import ins_interface

class ins_main(ins_interface):

    def __init__(self, ins_vars, **kwargs):

        """
        Constructor for the ins unit

        Arguments
        ---------

        ins_vars : list
                List of string for field variables required by ins unit
               
                ins_vars[0] --> velocity
                ins_vars[1] --> RHS from the previous time step
                ins_vars[2] --> divergence
                ins_vars[3] --> pressure


        **kwargs : Dictionary of keyword arguments

        'time_stepping' keyword refers to the time advancement scheme to be used

        kwargs['time_stepping'] = 'ab2' --> default
                                = 'euler'

        """

        from flowx.ins.solvers.projection import predictor, predictor_AB2, predictor_RK3, corrector, divergence
        from flowx.ins.solvers.stats import stats
        from flowx.ins.solvers.mass_balance import get_qin, get_qout, rescale_velocity, get_convvel, update_outflow_bc

        self.velc = ins_vars[0]
        self.hvar = ins_vars[1]
        self.divv = ins_vars[2]
        self.pres = ins_vars[3]

        self._time_stepping = 'ab2'

        if 'time_stepping' in kwargs: self._time_stepping = kwargs['time_stepping']

        if self._time_stepping is 'euler':
            self.predictor = predictor
        elif self._time_stepping is 'ab2':
            self.predictor = predictor_AB2
  
        self.divergence = divergence
        self.corrector = corrector
        self.stats = stats

        self.get_qin = get_qin
        self.get_qout = get_qout
        self.rescale_velocity = rescale_velocity
        self.get_convvel = get_convvel
        self.update_outflow_bc = update_outflow_bc

        return

    def advance(self, poisson, imbound, gridc, gridx, gridy, scalars, particles):

        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations
 
        Arguments
        ---------
        poisson : object
            Object for the poisson solver

        imbound : object
            Object for the immersed boundary unit

        gridc : object
          Grid object for cell centered variables

        gridx : object
          Grid object for x-face variables

        gridy : object
          Grid object for y-face variables

        scalars: object
           Scalars object to access time-step and Reynold number

        particles: object
           Object containing immersed boundary information
        """

        # Compute mass in
        Qin =  self.get_qin(gridx, self.velc) + self.get_qin(gridy, self.velc)

        # Update BC for predictor step
        self.update_outflow_bc(gridx, self.velc, scalars.variable['dt'])
        self.update_outflow_bc(gridy, self.velc, scalars.variable['dt'])

        # Calculate predicted velocity: u* = dt*H(u^n)       
        self.predictor(gridx, gridy, self.velc, self.hvar, scalars.variable['Re'], scalars.variable['dt'])
 
        # Immersed boundary forcing
        imbound.force_flow(gridx, gridy, scalars, particles)

        gridx.fill_guard_cells(self.velc)
        gridy.fill_guard_cells(self.velc)

        # Calculate RHS for the pressure Poission solver div(u)/dt
        self.divergence(gridc, gridx, gridy, self.velc, self.divv, ifac = scalars.variable['dt'])
        gridc.fill_guard_cells(self.divv)

        # Compute mass out
        Qout =  self.get_qout(gridx, self.velc) + self.get_qout(gridy, self.velc)

        # Rescale velocity to balance mass
        self.rescale_velocity(gridx, self.velc, Qin, Qout)
        self.rescale_velocity(gridy, self.velc, Qin, Qout)

        # Update BC for corrector step
        self.update_outflow_bc(gridx, self.velc, scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])
        self.update_outflow_bc(gridy, self.velc, scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])

        # Solve pressure Poisson equation
        scalars.stats['ites'], scalars.stats['res'] = poisson.solve_poisson(gridc)

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        self.corrector(gridc, gridx, gridy, self.velc, self.pres, scalars.variable['dt'])
        gridx.fill_guard_cells(self.velc)
        gridy.fill_guard_cells(self.velc)
   
        # Calculate divergence of the corrected velocity to display stats
        self.divergence(gridc, gridx, gridy, self.velc, self.divv)
        gridc.fill_guard_cells(self.divv)

        # Calculate stats
        scalars.stats.update(self.stats(gridc, gridx, gridy, self.velc, self.pres, self.divv))

        return
