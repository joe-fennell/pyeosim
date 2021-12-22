import os
import numpy as np
import json


class GenericTransformer(object):
    """A generic base class.

    Child subclasses should implement a '_set_steps' method that assigns a list
    of tuples to the 'steps' attribute.
    """

    def __init__(self):
        self.store_steps = False
        self.step_outputs = {}
        self.steps = []
        self._fitted = False  # flag to check if fixed pattern noise generated

    def fit(self, signal):
        """Sets all steps in the model.

        Args:
            signal: An xarray signal array
        """
        self._set_steps()
        self._fitted = True

    def transform(self, signal):
        """Applies all steps specificied in the 'steps' list.

        Args:
            signal: An xarray signal array

        Returns:
            new_signal: transformed values
        """
        try:
            signal = signal.transpose('y', 'x', ...)
        except ValueError:
            raise ValueError('x or y coordinate missing')
        if self._fitted:
            return self._apply_steps(signal)
        else:
            raise RuntimeError('Transformer Instance not fitted')

    def fit_transform(self, signal):
        """Calls 'fit' then 'transform' method.

        Args:
            signal: An xarray signal array

        Returns:
            new_signal: transformed values
        """
        self.fit(signal)
        return self.transform(signal)

    def get_steps(self):
        """
        Returns list of all steps in transformer
        """
        return [x[0] for x in self.steps]

    def get_params(self, numeric_only=False):
        """Returns the configuration params for the transfomer.

        Args:
            numeric_only (bool): if True returns only numeric type parameters

        Returns:
            A dictionary of parameter values keyed by name
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

    def apply_step(self, signal, step_name):
        """Apply a named step to a signal array

        Args:
            signal: An xarray signal array
            step_name: name of step to apply

        Returns:
            new_signal: transformed values
        """
        step = self.steps[step_name]
        return signal.pipe(step[1], **step[2])

    def _apply_steps(self, signal):
        # metadata copying and storag
        meta = signal.attrs.copy()
        _signal = signal.copy()
        if self.store_steps:
            for step in self.steps:
                _signal = _signal.pipe(step[1], **step[2]).compute()
                self.step_outputs[step[0]] = _signal.copy()
        else:
            for step in self.steps:
                _signal = _signal.pipe(step[1], **step[2])
            _signal = _signal.compute()

        name = self.__class__.__name__
        k = '{}_sensor_simulation'.format(name)
        meta[k] = json.dumps(self.get_params(),
                             default=lambda o: '<Not stored>')
        # meta[k] = str(self.get_params())
        _signal.attrs = meta
        return _signal

    def get_steps_as_latex(self, filepath=None):
        """Transforms steps to a latex table.

        Args:
            filepath (:obj:`str`, optional): destination to write latex file.

        Returns:
            Returns a LaTex string if filepath is not set
        """
        out = ''

        def add_line(x):
            line = x+'\n'
            line = line.replace('_', '\_')
            out += line
            # with open(filename, 'a') as f:
            #     f.write(line)

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
                # array-like
                val = list(val.values)
                val = '[{:.2f}, ...]'.format(np.ravel(val)[0])
            except AttributeError:
                # float, int or numeric string
                try:
                    val = float(val)
                    if (val > 100) | (val < .01):
                        val = '{:.2e}'.format(val)
                    else:
                        val = '{:.2f}'.format(val)
                # finally, just replace with NA
                except:
                    val = 'N/A'

            return '{} = '.format(name) + val

        def fmt_empty_step(val):
            return ' & & & ' + val + '\\\\'


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
        if filepath is not None:
            if os.path.exists(filepath):
                os.remove(filepath)
            with open(filename, 'a') as f:
                f.write(out)
        else:
            return out
