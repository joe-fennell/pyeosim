class GenericTransformer(object):
    """
    A generic base class. Child subclasses should implement a 'set_steps'
    method that assigns a list of tuples to the 'named_steps' attribute.
    """

    def __init__(self):
        self.store_steps = False
        self.step_outputs = {}

    def fit(self, signal):
        pass

    def transform(self, signal):
        meta = signal.attrs.copy()
        new = self._apply_steps(signal)
        name = self.__class__.__name__
        k = '{}_sensor_simulation'.format(name)
        meta[k] = str(self.get_params())
        new.attrs = meta
        return new

    def get_steps(self):
        """
        Returns list of all steps in transformer
        """
        return [x[0] for x in self.steps]

    def get_params(self):
        """
        Returns the configuration params for the transfomer
        """
        params = self.__dict__.copy()
        for k in ['steps', 'step_outputs', 'store_steps']:
            params.pop(k)
        return params

    def _apply_steps(self, signal):
        _signal = signal.copy()
        if self.store_steps:
            for step in self.steps:
                _signal = _signal.pipe(step[1], *step[2]).compute()
                self.step_outputs[step[0]] = _signal
            return _signal
        else:
            for step in self.steps:
                _signal = _signal.pipe(step[1], *step[2])
            return _signal.compute()
