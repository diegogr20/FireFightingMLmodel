class Particle:
    def __init__(self, **kwargs):
        self._properties = {"state","x", "tmem", "ux", "factor"}
        self._props = {"state", "x", "ux"}
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

