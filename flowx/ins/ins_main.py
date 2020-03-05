"""Module for incompressible Navier Stokes equations"""

from flowx.ins.ins_interface import ins_interface

class ins_main(ins_interface):

    def __init__(self, ins_vars=None, ins_info=None):

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


        ins_info : Dictionary of keyword arguments

        'time_stepping' keyword refers to the time advancement scheme to be used

        ins_info['time_stepping'] = 'ab2' --> default
                                  = 'euler'

        """

        from flowx.ins.solvers.projection import predictor, predictor_AB2, predictor_RK3, corrector, divergence
        from flowx.ins.solvers.stats import stats
        from flowx.ins.solvers.mass_balance import get_qin, get_qout, rescale_velocity, get_convvel, update_outflow_bc

        self._velc = 'stub'
        self._hvar = 'stub'
        self._divv = 'stub'
        self._pres = 'stub'

        self._time_stepping = 'ab2'

        if ins_info and 'time_stepping' in ins_info: self._time_stepping = ins_info['time_stepping']

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

        if ins_vars:
            self._velc = ins_vars[0]
            self._hvar = ins_vars[1]
            self._divv = ins_vars[2]
            self._pres = ins_vars[3]

        else:
            print('Warning: Incomp NS unit is a stub, any call to its methods will result in an error.') 

        return

    def advance(self, poisson, imbound, domain_data_struct):

        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations
 
        Arguments
        ---------
        poisson : object
            Object for the poisson solver

        imbound : object
            Object for the immersed boundary unit

        domain_data_struct : object list
           [gridc, gridx, gridy, scalars, particles]
        """

        _gridc = domain_data_struct[0]
        _gridx = domain_data_struct[1]
        _gridy = domain_data_struct[2]
        _scalars = domain_data_struct[3]

        # Compute mass in
        _Qin =  self._get_qin(_gridx, self._velc) + self._get_qin(_gridy, self._velc)

        # Update BC for predictor step
        self._update_outflow_bc(_gridx, self._velc, _scalars.variable['dt'])
        self._update_outflow_bc(_gridy, self._velc, _scalars.variable['dt'])

        # Calculate predicted velocity: u* = dt*H(u^n)       
        self._predictor(_gridx, _gridy, self._velc, self._hvar, _scalars.variable['Re'], _scalars.variable['dt'])
 
        # Immersed boundary forcing
        imbound.force_flow(domain_data_struct)

        _gridx.fill_guard_cells(self._velc)
        _gridy.fill_guard_cells(self._velc)

        # Calculate RHS for the pressure Poission solver div(u)/dt
        self._divergence(_gridc, _gridx, _gridy, self._velc, self._divv, ifac = _scalars.variable['dt'])
        _gridc.fill_guard_cells(self._divv)

        # Compute mass out
        _Qout =  self._get_qout(_gridx, self._velc) + self._get_qout(_gridy, self._velc)

        # Rescale velocity to balance mass
        self._rescale_velocity(_gridx, self._velc, _Qin, _Qout)
        self._rescale_velocity(_gridy, self._velc, _Qin, _Qout)

        # Update BC for corrector step
        self._update_outflow_bc(_gridx, self._velc, _scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])
        self._update_outflow_bc(_gridy, self._velc, _scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])

        # Solve pressure Poisson equation
        _scalars.stats['ites'], _scalars.stats['res'] = poisson.solve_poisson(_gridc)

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        self._corrector(_gridc, _gridx, _gridy, self._velc, self._pres, _scalars.variable['dt'])
        _gridx.fill_guard_cells(self._velc)
        _gridy.fill_guard_cells(self._velc)
   
        # Calculate divergence of the corrected velocity to display stats
        self._divergence(_gridc, _gridx, _gridy, self._velc, self._divv)
        _gridc.fill_guard_cells(self._divv)

        # Calculate stats
        _scalars.stats.update(self._stats(_gridc, _gridx, _gridy, self._velc, self._pres, self._divv))

        return
