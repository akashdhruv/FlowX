"""Module for incompressible Navier Stokes equations"""

from flowx.ins.ins_interface import ins_interface

class ins_main(ins_interface):

    def __init__(self, poisson=None, imbound=None, domain_data_struct=[None]*5, ins_vars=[None]*4, ins_info=None):

        """
        Constructor for the ins unit

        Arguments
        ---------

        poisson : object
            Object for the poisson solver

        imbound : object
            Object for the immersed boundary unit

        domain_data_struct : object list
           [gridc, gridx, gridy, scalars, particles]

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


        [self._gridc, self._gridx, self._gridy, self._scalars, self._particles] = domain_data_struct
        [self._velc, self._hvar, self._divv, self._pres] = ins_vars

        self._options = {'time_stepping' : 'ab2'}

        if ins_info:
            for key in ins_info: self._options[key] = ins_info[key]


        self._predictor_type = {'euler' : predictor, 'ab2': predictor_AB2}
        self._predictor = self._predictor_type[self._options['time_stepping']]

        self._divergence = divergence
        self._corrector = corrector
        self._stats = stats

        self._get_qin = get_qin
        self._get_qout = get_qout
        self._rescale_velocity = rescale_velocity
        self._get_convvel = get_convvel
        self._update_outflow_bc = update_outflow_bc

        self._imbound = imbound
        self._poisson = poisson

        self._ins_advance = self._advance_

        if None in domain_data_struct or None in ins_vars:
            self._ins_advance = self._advance_stub
            print('Warning: Incomp NS unit is a stub.') 

    def advance(self):
        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations 
        """

        self._ins_advance()
  
        return

    def _advance_stub(self):
        """
        Stub subroutine for the fractional step explicit time advancement of Navier Stokes equations 
        """

        return

    def _advance_(self):

        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations 
        """

        # Compute mass in
        _Qin =  self._get_qin(self._gridx, self._velc) + self._get_qin(self._gridy, self._velc)

        # Update BC for predictor step
        self._update_outflow_bc(self._gridx, self._velc, self._scalars.variable['dt'])
        self._update_outflow_bc(self._gridy, self._velc, self._scalars.variable['dt'])

        # Calculate predicted velocity: u* = dt*H(u^n)       
        self._predictor(self._gridx, self._gridy, self._velc, self._hvar, \
                        self._scalars.variable['Re'],self._scalars.variable['dt'])
 
        # Immersed boundary forcing
        self._imbound.force_flow()

        self._gridx.fill_guard_cells(self._velc)
        self._gridy.fill_guard_cells(self._velc)

        # Calculate RHS for the pressure Poission solver div(u)/dt
        self._divergence(self._gridc, self._gridx, self._gridy, self._velc, self._divv, ifac = self._scalars.variable['dt'])
        self._gridc.fill_guard_cells(self._divv)

        # Compute mass out
        _Qout =  self._get_qout(self._gridx, self._velc) + self._get_qout(self._gridy, self._velc)

        # Rescale velocity to balance mass
        self._rescale_velocity(self._gridx, self._velc, _Qin, _Qout)
        self._rescale_velocity(self._gridy, self._velc, _Qin, _Qout)

        # Update BC for corrector step
        self._update_outflow_bc(self._gridx, self._velc, self._scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])
        self._update_outflow_bc(self._gridy, self._velc, self._scalars.variable['dt'], convvel=[0.0,0.0,0.0,0.0])

        # Solve pressure Poisson equation
        self._scalars.stats['ites'], self._scalars.stats['res'] = self._poisson.solve_poisson()

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        self._corrector(self._gridc, self._gridx, self._gridy, self._velc, self._pres, self._scalars.variable['dt'])
        self._gridx.fill_guard_cells(self._velc)
        self._gridy.fill_guard_cells(self._velc)
   
        # Calculate divergence of the corrected velocity to display stats
        self._divergence(self._gridc, self._gridx, self._gridy, self._velc, self._divv)
        self._gridc.fill_guard_cells(self._divv)

        # Calculate stats
        self._scalars.stats.update(self._stats(self._gridc, self._gridx, self._gridy, self._velc, self._pres, self._divv))

        return
