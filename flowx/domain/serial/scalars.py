class Scalars(object):
    """
    Class to store and advance scalar data such as simulation time, 
    time-step, Reynolds number, etc.
    """

    def __init__(self,**kwargs):
        """
        This constuctor allows user to create scalar variables 
        specific to their simulation. Essential variables are defined 
        by default. The class stores two dictionaries,

        variable - to store values for time-step, Reynolds number, tmax, etc.

        stats - to store statistics to monitor simulation progress. 
        """

        self._set_default_values()
        self._set_user_values(**kwargs)

    def _set_default_values(self):
        """
        Private subroutine to set default values
        """

        var_list = ['to', 'tmax', 'time', 'dt', 'nstep']
        val_list = [0.0, 0.0, 0.0, 1.0, 0]

        self.variable = dict(zip(var_list, val_list))
        self.stats = dict()

    def _set_user_values(self,**kwargs):
        """
        Private subroutine to set user defined values
        """
        for key,value in kwargs.items():
            self.variable[key] = value

    def advance(self):
        """
        Subroutine to advance the simulation time and time-step
        """
        self.variable['time'] += self.variable['dt']
        self.variable['nstep'] += 1
