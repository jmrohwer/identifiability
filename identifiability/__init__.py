#!/usr/bin/env python

"""
Custom functions for identifiability analysis to calculate 
and plot confidence intervals based on a profile-likelihood analysis. Adapted 
from lmfit, with custom functions to select the range for parameter scanning and 
for plotting the profile likelihood.
"""

from collections import OrderedDict
from lmfit.minimizer import MinimizerException, Minimizer, MinimizerResult
from lmfit.model import ModelResult
import numpy as np
import scipy as sp
import math
from matplotlib import pyplot as plt
try:
    from multiprocessing_on_dill import Pool
except ModuleNotFoundError:
    print(
        """Info: module 'multiprocessing_on_dill' is not installed.
Parameter estimation of kinetic models with PySCeS using the CVODE solver won't work!
"""
    )
    from multiprocessing import Pool

from .version import __version__

CONF_ERR_GEN = 'Cannot determine Confidence Intervals'
CONF_ERR_NVARS = '%s with < 2 variables' % CONF_ERR_GEN

def _copy_vals(params):
    """Save values/stderrs of parameters in a temporary dictionary."""
    tmp_params = {}
    for para_key in params:
        tmp_params[para_key] = (params[para_key].value,
                                params[para_key].stderr)
    return tmp_params

class ConfidenceInterval:
    """Class used to calculate the confidence interval."""

    def __init__(self, minimizer, result, p_names=None):
        """Initialize the ConfidenceInterval class.

        Parameters
        ----------
        minimizer : Minimizer
            The minimizer to use, holding objective function.
        result : MinimizerResult
            The result of running minimize().
        p_names : list, optional
            Names of the parameters for which the CI is calculated. If None
            (default), the CI is calculated for every parameter.

        Raises
        ------
        MinimizerException
            If there are less than two variables.
        """

        assert isinstance(minimizer, Minimizer) or isinstance(minimizer, ModelResult), (
            'minimizer must be instance of `lmfit.minimizer.Minimizer` or `lmfit.model.ModelResult`'
        )
        assert isinstance(result, MinimizerResult) or isinstance(result, ModelResult), (
            'result must be instance of `lmfit.minimizer.MinimizerResult` or `lmfit.model.ModelResult`'
        )

        self.minimizer = minimizer
        self.result = result
        self.params = result.params.copy()
        self.org = _copy_vals(self.params)
        self.best_chi = result.chisqr

        if not p_names:
            p_names = [i for i in self.params if self.params[i].vary]
        self.p_names = p_names
        self.fit_params = [self.params[p] for p in self.p_names]

        self._traces_calculated = False
        self._k = 2   # degree of smoothing spline

        # check that there are at least 2 true variables!
        nvars = len([p for p in self.params.values() if p.vary])
        if nvars < 2:
            raise MinimizerException(CONF_ERR_NVARS)

        self.trace_dict = {i: {} for i in self.p_names}

    def calc_all_ci(
        self,
        limits=0.5,
        points=11,
        prob=0.95,
        method='leastsq',
        log=False,
        recalc=False,
        mp=True,
    ):
        """Calculate all confidence intervals.

        Parameters
        ----------
        limits : float, optional
            The limits (as a fraction of the original parameter value) within which
            to vary the parameters for identifiability analysis (default is 0.5).
            If ``log=False``, the parameter is varied from p*limits to p*(2 - limits),
            where p is the original value.
            If ``log=True``, the parameter is varied from p*limits to p/limits.
        points : int, optional
            The number of points for which to calculate the profile likelihood over
            the given parameter range.
        prob : float, optional
            The probability for the confidence interval (<1). If None,
            the default is 0.95 (95 % confidence interval).
        method : str, optional
            The lmfit mimimize() method to use (default='leastsq')
        log : bool, optional
            Whether to vary the parameter in a log (True) or a linear (False,
            default) scale.
        recalc : bool, optional
            Whether to recalculate the traces and splines (default=False). If False,
            the existing traces are used (useful for calculating the CIs for a different
            probability). If True, the traces are recalculated. This is useful
            if the `limits`, `points`, `method` or `log`  parameters are changed.
        mp : bool, optional
            Run the optimization in parallel using ``multiprocessing`` (default=True)

        Returns
        -------
        output : dict
            A dictionary containing a list of ``(lower, upper)``-tuples containing
            the confidence bounds for each parameter.
        """

        assert (
            (type(prob) == float) & (prob > 0) & (prob < 1)
        ), 'Please provide a probability value between 0 and 1.'
        self.prob = prob
        self.method = method
        self.log = log
        self.ci_values = OrderedDict()
        self.threshold = self._calc_threshold()

        if not self._traces_calculated or recalc:
            self._populate_traces(limits, points, mp)

        for p in self.p_names:
            self.ci_values[p] = self._process_ci(p)

        return self.ci_values

    def _populate_traces(self, limits, points, mp):
        if mp:
            proc_pool = Pool()
            arl = []

        results = []

        for para in self.p_names:
            if isinstance(para, str):
                para = self.params[para]

            if self.log:
                para_vals = np.logspace(
                    np.log10(para.value * limits), np.log10(para.value / limits), points,
                )
            else:
                para_vals = np.linspace(limits * para.value, (2 - limits) * para.value, points)

            para.vary = False
            self.trace_dict[para.name]['value'] = []
            self.trace_dict[para.name]['dchi'] = []
            self.trace_dict[para.name]['results'] = []

            for val in para_vals:
                self.trace_dict[para.name]['value'].append(val)
                if mp:
                    arl.append(proc_pool.apply_async(self._calc_dchi, args=(self, para, val)))
                else:
                    results.append(self.calc_dchi(para, val))

            para.vary = True
            self._reset_vals()

        if mp:
            arl[-1].wait()
            for ar in arl:
                results.append(ar.get())
            proc_pool.terminate()
            del proc_pool

        for (para, dchi, opt_res) in results:
            self.trace_dict[para.name]['dchi'].append(dchi)
            self.trace_dict[para.name]['results'].append(opt_res)
        self._traces_calculated = True


    def _process_ci(self, p_name):
        xx = self.trace_dict[p_name]['value']
        yy = self.trace_dict[p_name]['dchi']
        t = self.threshold
        spl = sp.interpolate.UnivariateSpline(xx, yy, k=self._k, s=0)
        if self.log:
            allx = np.logspace(np.log10(xx[0]), np.log10(xx[-1]), 20000)
        else:
            allx = np.linspace(xx[0], xx[-1], 20000)

        lo = allx[spl(allx) <= t][0]
        hi = allx[spl(allx) <= t][-1]

        # catch non-identifiable cases
        if lo == xx[0]:
            lo = np.nan
        if hi == xx[-1]:
            hi = np.nan
        return lo, hi

    def _reset_vals(self):
        """Reset parameter values to best-fit values."""
        for para_key in self.params:
            (self.params[para_key].value, self.params[para_key].stderr,) = self.org[
                para_key
            ]

    @staticmethod
    def _calc_dchi(ci_instance, para, val):
        """
        Static method to calculate the normalised delta chi-squared
        using multiprocessing.
        """
        para.vary = False
        para.value = val
        save_para = ci_instance.params[para.name]
        ci_instance.params[para.name] = para
        ci_instance.minimizer.prepare_fit(ci_instance.params)
        out = ci_instance.minimizer.minimize(method=ci_instance.method)
        dchi = ci_instance._dchi(ci_instance.result, out)
        ci_instance.params[para.name] = save_para
        para.vary = True
        return para, dchi, out

    def calc_dchi(self, para, val, restore=False):
        """
        Calculate the normalised delta chi-squared for 
        a given parameter value.
        """
        if restore:
            self._reset_vals()
        para.value = val
        save_para = self.params[para.name]
        self.params[para.name] = para
        self.minimizer.prepare_fit(self.params)
        out = self.minimizer.minimize(method=self.method)
        dchi = self._dchi(self.result, out)
        self.params[para.name] = save_para
        return para, dchi, out

    def _dchi(self, best_fit, new_fit):
        """
        Return the normalised delta chi-squared between the best fit
        and the new fit.
        """
        dchi = new_fit.chisqr / best_fit.chisqr - 1.0
        return dchi

    def _calc_threshold(self):
        """
        Return the threshold of the normalised chi-squared for 
        the given probability.
        """
        nfree = self.result.nfree
        nfix = 1
        threshold_scaled = sp.stats.chi2.ppf(self.prob, nfix)
        threshold = threshold_scaled * nfix / nfree
        return threshold

    def plot_ci(self, para, ax=None):
        """Plot profile likelihood with confidence interval for single parameter.

        Parameters
        ----------
        para : str
            The parameter name for which to plot the profile likelihood.
        ax: matplotlib.axes.Axes, optional
            Matplotlib Axes object to plot on. If None, a new figure and axes will be created.

        """

        assert para in self.p_names, 'para must be one of ' + str(self.p_names)
        if not ax:
            f, ax = plt.subplots()
        xx = self.trace_dict[para]['value']
        yy = self.trace_dict[para]['dchi']
        t = self.threshold
        spl = sp.interpolate.UnivariateSpline(xx, yy, k=self._k, s=0)
        allx = np.linspace(xx[0], xx[-1], 20000)
        ax.plot(xx, yy, '+')
        ax.plot(allx, spl(allx), '-', lw=1)
        ax.axhline(t, color='k', ls='--', lw=0.5)
        ax.axvline(self.params[para].value, color='k', ls='-', lw=0.5)
        lo, hi = self.ci_values[para]
        if np.isnan(lo):
            lo = ax.get_xlim()[0]
        if np.isnan(hi):
            hi = ax.get_xlim()[1]
        ax.axvspan(lo, hi, alpha=0.1, color='b')
        if self.log:
            ax.semilogx()
        ax.set_xlabel('Parameter value')
        ax.set_ylabel(r'$\chi^2\left/\chi^2_0\right. - 1$')
        ax.set_title(para)

    def plot_all_ci(self):
        """Plot profile likelihoods with confidence intervals for all parameters."""

        num = len(self.p_names)
        numcols = 3
        numrows = math.ceil(num / numcols)
        f, ax = plt.subplots(nrows=numrows, ncols=numcols, figsize=(9, 2.5 * numrows))
        for i in range(num):
            if num <= numcols:
                theax = ax[i]
            else:
                theax = ax[i // numcols, i % numcols]
            self.plot_ci(self.p_names[i], ax=theax)
        # remove empty axes
        if num % numcols != 0:
            empty = numcols - num % numcols
            for i in range(-empty, 0):
                if num <= numcols:
                    ax[i].set_visible(False)
                else:
                    ax[num // numcols, i].set_visible(False)
        f.tight_layout()

