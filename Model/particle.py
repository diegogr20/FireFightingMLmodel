class Particle:
    """Object to model a convection or a radiation particle.

    Params:
    ------
    state  (float): For each particle, a value 0 or 1 ("cold" or "ignited").
                    If "ignited", it moves until it dies according to a decay.
    tmem   (float): Timescale of the decay (s) for each particle - for now, this
                    reflects the QMAX
                    
    x      (float): x-position of each particle
    y      (float): y-position of each particle

    factor (float): From 0 to 1, according to a Gaussian, centred at max burning
    ux     (float): mean wind speed in X-dir (later, it can be f(Z))
    vy     (float): mean wind speed in Y-dir


    Attributes:
    -----------
    props_name: all particle's properties name
                READ ONLY
    values: a dict containing all values for each key in props
            READ ONLY


    Methods:
    --------

    update(self, **kwargs): update the particle's value using kwargs dict

    get_from_keys(self, keys): return a tuple with the values of each property
                               in the iterable 'keys' IN THE GIVEN ORDER.


    time, with spread so as to account for IGNTIME & BURNDUR.
    In this version, all cells emit equal number of particles, but they
    could have different MEMORY to account for different QMAX
    """

    def __init__(self, **kwargs):
        self._properties = {"state","x", "y", "tmem", "ux", "vy", "factor"}
        self._props = {"state", "x", "y", "ux", "vy"}
        self._values = {props: 0 for props in self._properties}
        self.update(**kwargs)

    @property
    def values(self):
        """Get a key: value dictionnary containning the properties
        of the particle."""
        return self._values

    def update(self, **kwargs):
        for key, value in kwargs.items():
            if key in self._properties:
                self._values[key] = value

    def get_from_keys(self, keys):
        val = ()
        for key in keys:
            val += (self[key], )
        return val

    def __getattr__(self, key):
        if key in self._properties:
            return self._values[key]
        else:
            raise AttributeError(f"particle does not have property '{key}'")

    def __setattr__(self, key, value):
        if key in {"_properties", "_values", "_props"}:
            object.__setattr__(self, key, value)
        elif key in self._properties:
            self._values[key] = value

    def __getitem__(self, key):
        if key in self._properties:
            return self._values[key]
        raise KeyError(f"particle doesn't have property {key}.")

