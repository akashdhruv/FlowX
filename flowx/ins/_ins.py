"""Module for incompressible Navier Stokes equations"""

import time

from . import _interface


class IncompNS(object):
    def __init__(
        self,
        poisson=None,
        imbound=None,
        domain_data_struct=[None] * 5,
        ins_vars=[None] * 5,
        ins_info=None,
    ):
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
                ins_vars[4] --> delta pressure

        ins_info : Dictionary of keyword arguments

        'time_stepping' keyword refers to the time advancement scheme to be used

        ins_info['time_stepping'] = 'ab2' --> default
                                  = 'euler'
        """
        # ----------Create images for other unit objects and variables----------------------------
        (
            self._gridc,
            self._gridx,
            self._gridy,
            self._scalars,
            self._particles,
        ) = domain_data_struct
        self._velc, self._hvar, self._divv, self._pres, self._delp = ins_vars

        self._imbound = imbound
        self._poisson = poisson

        # ----------Setup default parameters for the current unit---------------------------------
        self._options = {"time_stepping": "ab2", "pressure_correct": True}

        self._predictor_type = {
            "euler": _interface.predictor_euler,
            "ab2": _interface.predictor_ab2,
        }

        # ----------------------Read user parameters----------------------------------------------
        if ins_info:
            for key in ins_info:
                self._options[key] = ins_info[key]

        # ---------------------Setup unit based on user/default parameters----------------------
        if (
            None in domain_data_struct
            or None in ins_vars
            or imbound is None
            or poisson is None
        ):
            self._ins_advance = self._advance_stub
            print(
                "Warning: Incomp NS unit is a stub because one or more parameters were not supplied."
            )

        elif (
            self._gridc.nblocks > 1
            or self._gridx.nblocks > 1
            or self._gridy.nblocks > 1
        ):
            self._ins_advance = self._advance_stub
            print("Warning: Incomp NS unit is a stub because nblocks > 1.")

        else:
            self._ipres = self._options["pressure_correct"]
            self._predictor = self._predictor_type[self._options["time_stepping"]]

            self._divergence = _interface.divergence
            self._corrector = _interface.corrector
            self._stats = _interface.stats

            self._get_qin = _interface.get_qin
            self._get_qout = _interface.get_qout
            self._rescale_velocity = _interface.rescale_velocity
            self._get_convvel = _interface.get_convvel
            self._update_outflow_bc = _interface.update_outflow_bc

            self._ins_advance = self._advance

            self._ustar_bc = self._gridx.bc_type[self._velc].copy()
            self._vstar_bc = self._gridy.bc_type[self._velc].copy()

            self._ucorr_bc = self._ustar_bc.copy()
            self._vcorr_bc = self._vstar_bc.copy()

            self._pres_bc = self._gridc.bc_type[self._pres].copy()

            (
                self._ucorr_bc[0],
                self._ucorr_bc[1],
                self._vcorr_bc[2],
                self._vcorr_bc[3],
            ) = [None] * 4

            if "dirichlet" in self._pres_bc:
                self._rescaleVelocity = _interface.rescale_velocity_stub

        return

    def advance(self):
        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations
        """
        self._ins_advance()

    def _advance_stub(self):
        """
        Stub subroutine for the fractional step explicit time advancement of Navier Stokes equations
        """
        pass

    def _advance(self):
        """
        Subroutine for the fractional step explicit time advancement of Navier Stokes equations
        """
        time_ins_start = time.time()

        # Update BC for predictor step
        self._gridx.update_bc_type({self._velc: self._ustar_bc})
        self._gridy.update_bc_type({self._velc: self._vstar_bc})

        # Compute mass in
        _Qin = self._get_qin(self._gridx, self._velc) + self._get_qin(
            self._gridy, self._velc
        )
        self._scalars.stats["qin"] = _Qin

        # Calculate outflow BC
        self._update_outflow_bc(self._gridx, self._velc, self._scalars.dt)
        self._update_outflow_bc(self._gridy, self._velc, self._scalars.dt)

        # Calculate predicted velocity: u* = dt*H(u^n)
        self._predictor(
            self._gridc,
            self._gridx,
            self._gridy,
            self._velc,
            self._hvar,
            self._pres,
            self._scalars.Re,
            self._scalars.dt,
            self._ipres,
        )

        # Immersed boundary forcing
        self._imbound.force_flow()

        self._gridx.fill_guard_cells(self._velc)
        self._gridy.fill_guard_cells(self._velc)

        # Compute mass out
        _Qout = self._get_qout(self._gridx, self._velc) + self._get_qout(
            self._gridy, self._velc
        )
        self._scalars.stats["qout"] = _Qout

        # Rescale velocity to balance mass
        self._rescale_velocity(self._gridx, self._velc, _Qin, _Qout)
        self._rescale_velocity(self._gridy, self._velc, _Qin, _Qout)

        # Calculate RHS for the pressure Poission solver div(u)/dt
        self._divergence(
            self._gridc,
            self._gridx,
            self._gridy,
            self._velc,
            self._divv,
            ifac=self._scalars.dt,
        )
        self._gridc.fill_guard_cells(self._divv)

        # Solve pressure Poisson equation
        time_poisson_begin = time.time()
        self._scalars.stats["ites"], self._scalars.stats["res"] = self._poisson.solve()
        time_poisson_end = time.time()
        self._scalars.stats["poisson_time"] = time_poisson_end - time_poisson_begin

        # Update BC for corrector step
        self._gridx.update_bc_type({self._velc: self._ucorr_bc})
        self._gridy.update_bc_type({self._velc: self._vcorr_bc})

        # Calculate corrected velocity u^n+1 = u* - dt * grad(P)
        self._corrector(
            self._gridc,
            self._gridx,
            self._gridy,
            self._velc,
            self._pres,
            self._delp,
            self._scalars.dt,
            self._ipres,
        )
        self._gridx.fill_guard_cells(self._velc)
        self._gridy.fill_guard_cells(self._velc)

        # Calculate divergence of the corrected velocity to display stats
        self._divergence(self._gridc, self._gridx, self._gridy, self._velc, self._divv)
        self._gridc.fill_guard_cells(self._divv)

        # Calculate total INS time
        time_ins_end = time.time()
        self._scalars.stats["ins_time"] = time_ins_end - time_ins_start

        # Calculate stats
        self._scalars.stats.update(
            self._stats(
                self._gridc,
                self._gridx,
                self._gridy,
                self._velc,
                self._pres,
                self._divv,
            )
        )
