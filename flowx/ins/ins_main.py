"""Module for incompressible Navier Stokes equations"""

from flowx.ins.ins_interface import ins_interface

class ins_main(ins_interface):

    def __init__(self, poisson=None, imbound=None, domain_data_struct=[None]*5, ins_vars=[None]*5, ins_info=None):

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
        from flowx.ins.solvers.mass_balance import get_qin, get_qout, rescale_velocity, get_convvel, update_outflow_bc, \
                                                   rescale_velocity_stub


        self._gridc, self._gridx, self._gridy, self._scalars, self._particles = domain_data_struct
        self._velc, self._hvar, self._divv, self._pres, self._delp = ins_vars

        self._options = {'time_stepping' : 'ab2', 'pressure_correct' : True}

        if ins_info:
            for key in ins_info: self._options[key] = ins_info[key]

        self._ipres = self._options['pressure_correct']

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

        if(self._velc is not None and self._pres is not None):

            self._ustar_bc = self._gridx.bc_type[self._velc].copy()
            self._vstar_bc = self._gridy.bc_type[self._velc].copy()
        
            self._ucorr_bc = self._ustar_bc.copy()
            self._vcorr_bc = self._vstar_bc.copy()

            self._pres_bc = self._gridc.bc_type[self._pres].copy() 

            self._ucorr_bc[0], self._ucorr_bc[1], self._vcorr_bc[2], self._vcorr_bc[3] = [None]*4
            
            if 'dirichlet' in self._pres_bc : self._rescale_velocity = rescale_velocity_stub

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

        import time

        time_ins_start = time.time()

        # Update BC for predictor step
        self._gridx.update_bc_type({self._velc: self._ustar_bc})
        self._gridy.update_bc_type({self._velc: self._vstar_bc})

        # Compute mass in
        _Qin =  self._get_qin(self._gridx, self._velc) + self._get_qin(self._gridy, self._velc)
        self._scalars.stats['qin'] = _Qin

        # Calculate outflow BC
        self._update_outflow_bc(self._gridx, self._velc, self._scalars.dt)
        self._update_outflow_bc(self._gridy, self._velc, self._scalars.dt)

        # Calculate predicted velocity: u* = dt*H(u^n)       
        self._predictor(self._gridc, self._gridx, self._gridy, self._velc, self._hvar, self._pres, \
                        self._scalars.Re, self._scalars.dt, self._ipres) 

        # Immersed boundary forcing
        self._imbound.force_flow()

        self._gridx.fill_guard_cells(self._velc)
        self._gridy.fill_guard_cells(self._velc)

        # Compute mass out
        _Qout =  self._get_qout(self._gridx, self._velc) + self._get_qout(self._gridy, self._velc)
        self._scalars.stats['qout'] = _Qout

        # Rescale velocity to balance mass
        self._rescale_velocity(self._gridx, self._velc, _Qin, _Qout)
        self._rescale_velocity(self._gridy, self._velc, _Qin, _Qout)

        # Calculate RHS for the pressure Poission solver div(u)/dt
        self._divergence(self._gridc, self._gridx, self._gridy, self._velc, self._divv, ifac = self._scalars.dt)
        self._gridc.fill_guard_cells(self._divv)

        # Solve pressure Poisson equation
        time_poisson_begin = time.time()
        self._scalars.stats['ites'], self._scalars.stats['res'] = self._poisson.solve_poisson()
        time_poisson_end = time.time()
        self._scalars.stats['poisson_time'] = time_poisson_end - time_poisson_begin

        # Update BC for corrector step
        self._gridx.update_bc_type({self._velc: self._ucorr_bc})
        self._gridy.update_bc_type({self._velc: self._vcorr_bc})

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P) 
        self._corrector(self._gridc, self._gridx, self._gridy, self._velc, self._pres, self._delp, \
                        self._scalars.dt, self._ipres)
        self._gridx.fill_guard_cells(self._velc)
        self._gridy.fill_guard_cells(self._velc)
   
        # Calculate divergence of the corrected velocity to display stats
        self._divergence(self._gridc, self._gridx, self._gridy, self._velc, self._divv)
        self._gridc.fill_guard_cells(self._divv)

        # Calculate total INS time
        time_ins_end = time.time()
        self._scalars.stats['ins_time'] = time_ins_end - time_ins_start

        # Calculate stats
        self._scalars.stats.update(self._stats(self._gridc, self._gridx, self._gridy, self._velc, self._pres, self._divv))

        return
