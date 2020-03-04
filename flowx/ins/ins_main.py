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

        self._velc = ins_vars[0]
        self._hvar = ins_vars[1]
        self._divv = ins_vars[2]
        self._pres = ins_vars[3]

        self._time_stepping = 'ab2'

        if 'time_stepping' in kwargs: self._time_stepping = kwargs['time_stepping']

        if self._time_stepping is 'euler':
            self._predictor = predictor
        elif self._time_stepping is 'ab2':
            self._predictor = predictor_AB2
  
        self._divergence = divergence
        self._corrector = corrector
        self._stats = stats

        self._get_qin = get_qin
        self._get_qout = get_qout
        self._rescale_velocity = rescale_velocity
        self._get_convvel = get_convvel
        self._update_outflow_bc = update_outflow_bc

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
        Qin =  self._get_qin(gridx, self._velc) + self._get_qin(gridy, self._velc)

        # Update BC for predictor step
        self._update_outflow_bc(gridx, self._velc, scalars.variable['dt'])
        self._update_outflow_bc(gridy, self._velc, scalars.variable['dt'])

        # Calculate predicted velocity: u* = dt*H(u^n)       
        self._predictor(gridx, gridy, self._velc, self._hvar, scalars.variable['Re'], scalars.variable['dt'])
 
        # Immersed boundary forcing
        imbound.force_flow(gridx, gridy, scalars, particles)

        gridx.fill_guard_cells(self._velc)
        gridy.fill_guard_cells(self._velc)

        # Calculate RHS for the pressure Poission solver div(u)/dt
        self._divergence(gridc, gridx, gridy, self._velc, self._divv, ifac = scalars.variable['dt'])
        gridc.fill_guard_cells(self._divv)

        # Compute mass out
        Qout =  self._get_qout(gridx, self._velc) + self._get_qout(gridy, self._velc)

        # Rescale velocity to balance mass
        self._rescale_velocity(gridx, self._velc, Qin, Qout)
        self._rescale_velocity(gridy, self._velc, Qin, Qout)

        # Update BC for corrector step
        self._update_outflow_bc(gridx, self._velc, scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])
        self._update_outflow_bc(gridy, self._velc, scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])

        # Solve pressure Poisson equation
        scalars.stats['ites'], scalars.stats['res'] = poisson.solve_poisson(gridc)

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        self._corrector(gridc, gridx, gridy, self._velc, self._pres, scalars.variable['dt'])
        gridx.fill_guard_cells(self._velc)
        gridy.fill_guard_cells(self._velc)
   
        # Calculate divergence of the corrected velocity to display stats
        self._divergence(gridc, gridx, gridy, self._velc, self._divv)
        gridc.fill_guard_cells(self._divv)

        # Calculate stats
        scalars.stats.update(self._stats(gridc, gridx, gridy, self._velc, self._pres, self._divv))

        return
