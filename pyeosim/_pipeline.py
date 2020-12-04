import os


class GenericTransformer(object):
    """
    A generic base class. Child subclasses should implement a 'set_steps'
    method that assigns a list of tuples to the 'named_steps' attribute.
    """

    def __init__(self):
        self.store_steps = False
        self.step_outputs = {}
        self._fitted = False  # flag to check if fixed pattern noise generated

    def fit(self, signal):
        self.set_steps()
        self._fitted = True
        pass

    def transform(self, signal):
        if self._fitted:
            return self._apply_steps(signal)
        else:
            raise RuntimeError('Transformer Instance not fitted')

    def fit_transform(self, signal):
        self.fit(signal)
        return self.transform(signal)

    def get_steps(self):
        """
        Returns list of all steps in transformer
        """
        return [x[0] for x in self.steps]

    def get_params(self, numeric_only=False):
        """
        Returns the configuration params for the transfomer
        """
        def only_numeric(params):
            def isnumeric(x):
                return (type(x) == float) | (type(x) == int)

            new = {}
            for k, v in params.items():
                if isnumeric(v):
                    new[k] = v
            return new

        params = self.__dict__.copy()
        for k in ['steps', 'step_outputs', 'store_steps']:
            params.pop(k)

        if numeric_only:
            return only_numeric(params)
        return params

    def _apply_steps(self, signal):
        # metadata copying and storag
        meta = signal.attrs.copy()
        _signal = signal.copy()
        if self.store_steps:
            for step in self.steps:
                _signal = _signal.pipe(step[1], **step[2]).compute()
                self.step_outputs[step[0]] = _signal
        else:
            for step in self.steps:
                _signal = _signal.pipe(step[1], **step[2])
            _signal = _signal.compute()

        name = self.__class__.__name__
        k = '{}_sensor_simulation'.format(name)
        meta[k] = str(self.get_params())
        _signal.attrs = meta
        return _signal

    def apply_step(self, signal, step_name):
        """
        Apply a named step to a signal array

        Parameters
        ----------
        signal : xarray.DataArray
            signal array
        step_name : str
            name of step to apply
        """
        step = self.steps[step_name]
        return signal.pipe(step[1], **step[2])

    def steps_to_latex(self, filename):
        """
        Transforms steps to a latex table.

        Parameters
        ----------
        filename : str
            path to write to
        """
        def add_line(x):
            line = x+'\n'
            line = line.replace('_', '\_')
            with open(filename, 'a') as f:
                f.write(line)

        def fmt_obj(obj):

            def parse_func(func):
                module = func.__module__
                name = func.__name__
                return module + '.' + name

            def parse_method(method):
                return str(method.__self__).split(
                    ' ')[0].replace('<', '')
            try:
                obj.__func__
                return parse_method(obj)

            except AttributeError:
                return parse_func(obj)

        def fmt_step(step, c, param_entry='None'):
            _str = '{} &'.format(c)
            _str += ' {} & '.format(step[0])
            _str += '{} & '.format(fmt_obj(step[1]))
            _str += '{}\\\\'.format(param_entry)
            return _str

        def fmt_params(name, val):
            try:
                val = list(val.values)
                val = '[{:.2f},...]'.format(val[0])
            except AttributeError:
                if type(val) != str:
                    if (val > 100) | (val < .01):
                        val = '{:.2e}'.format(val)
                    else:
                        val = '{:.2f}'.format(val)
            return '{} = '.format(name) + val

        def fmt_empty_step(val):
            return ' & & & ' + val + '\\\\'

        if os.path.exists(filename):
            os.remove(filename)

        # add all preamble
        add_line('\\begin{tabular}{llll}')
        add_line('\\toprule')
        add_line(' & \\bf{Stage} & \\bf{Function} & \\bf{Parameters} \\\\')
        add_line('\\midrule')
        # iterate all steps and write to a file
        c = 1
        for s in self.steps:
            n_params = len(s[2])
            if n_params == 0:
                add_line(fmt_step(s, c))
            if n_params >= 1:
                _pnames = list(s[2].keys())
                _pvals = list(s[2].values())
                _p = fmt_params(_pnames[0],
                                _pvals[0])
                add_line(fmt_step(s, c, _p))

            if n_params > 1:
                for i in range(1, n_params):
                    _p = fmt_params(_pnames[i], _pvals[i])
                    add_line(fmt_empty_step(_p))
            add_line('\\midrule')
            c += 1
        add_line('\\bottomrule')
        add_line('\\end{tabular}')
