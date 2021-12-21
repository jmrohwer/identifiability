#!/usr/bin/env python

"""
Custom functions for identifiability analysis to calculate 
and plot confidence intervals based on a profile-likelihood analysis. Adapted 
from lmfit, with custom functions to select the range for parameter scanning and 
for plotting the profile likelihood.
"""

__version__ = 0.1

from collections import OrderedDict
from lmfit.minimizer import MinimizerException
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

CONF_ERR_GEN = 'Cannot determine Confidence Intervals'
CONF_ERR_STDERR = '%s without sensible uncertainty estimates' % CONF_ERR_GEN
CONF_ERR_NVARS = '%s with < 2 variables' % CONF_ERR_GEN


class ConfidenceInterval:
    """Class used to calculate the confidence interval."""

    def __init__(self, minimizer, result, p_names=None, prob=None, log=False):
        self.minimizer = minimizer
        self.result = result
        self.params = result.params.copy()
        self.org = {}
        for para_key in self.params:
            self.org[para_key] = (
                self.params[para_key].value,
                self.params[para_key].stderr,
            )

        self.best_chi = result.chisqr

        if p_names is None:
            p_names = [i for i in self.params if self.params[i].vary]

        self.p_names = p_names
        self.fit_params = [self.params[p] for p in self.p_names]

        self.prob = prob
        self.log = log

        # check that there are at least 2 true variables!
        # check that all stderrs are sensible (including not None or NaN)

        for par in self.fit_params:
            if par.vary and (par.stderr is None or par.stderr is np.nan):
                raise MinimizerException(CONF_ERR_STDERR)
        nvars = len([p for p in self.params.values() if p.vary])
        if nvars < 2:
            raise MinimizerException(CONF_ERR_NVARS)

        self.trace_dict = {i: {} for i in self.p_names}

    def calc_all_ci(self, limits=0.5, points=11):
        """Calculate all confidence intervals."""
        self.ci_values = OrderedDict()

        for p in self.p_names:
            self.ci_values[p] = self.calc_ci(p, limits, points)

        return self.ci_values

    def calc_ci(self, para, limits, points):
        """Calculate the CI for a single parameter."""
        if isinstance(para, str):
            para = self.params[para]

        if self.log:
            para_vals = np.logspace(
                np.log10(para.value * limits), np.log10(para.value / limits), points,
            )
        else:
            para_vals = np.linspace(
                (1 - limits) * para.value, (1 + limits) * para.value, points
            )

        para.vary = False
        threshold = self.calc_threshold()
        self.trace_dict[para.name]['value'] = []
        self.trace_dict[para.name]['dchi'] = []
        self.trace_dict[para.name]['threshold'] = threshold

        for val in para_vals:
            self.trace_dict[para.name]['value'].append(val)
            self.trace_dict[para.name]['dchi'].append(self.calc_dchi(para, val))

        para.vary = True
        self.reset_vals()

        xx = self.trace_dict[para.name]['value']
        yy = self.trace_dict[para.name]['dchi']
        t = self.trace_dict[para.name]['threshold']
        spl = sp.interpolate.UnivariateSpline(xx, yy, k=2, s=0)
        if self.log:
            allx = np.logspace(np.log10(xx[0]), np.log10(xx[-1]), 20000)
        else:
            allx = np.linspace(xx[0], xx[-1], 20000)

        return (allx[spl(allx) <= t][0], allx[spl(allx) <= t][-1])

    def reset_vals(self):
        """Reset parameter values to best-fit values."""
        for para_key in self.params:
            (self.params[para_key].value, self.params[para_key].stderr,) = self.org[
                para_key
            ]

    def calc_dchi(self, para, val, restore=False):
        """
        Calculate the normalised delta chi-squared for 
        a given parameter value.
        """
        if restore:
            self.reset_vals()
        para.value = val
        save_para = self.params[para.name]
        self.params[para.name] = para
        self.minimizer.prepare_fit(self.params)
        out = self.minimizer.leastsq()
        dchi = self.dchi(self.result, out)
        self.params[para.name] = save_para
        return dchi

    def dchi(self, best_fit, new_fit):
        """
        Return the normalised delta chi-squared between the best fit
        and the new fit.
        """
        nfree = best_fit.nfree
        nfix = best_fit.nvarys - new_fit.nvarys
        dchi = new_fit.chisqr / best_fit.chisqr - 1.0
        return dchi

    def calc_threshold(self):
        """
        Return the threshold of the normalised chi-squared for 
        the given probability.
        """
        nfree = self.result.nfree
        nfix = 1
        threshold_scaled = sp.stats.chi2.ppf(self.prob, nfix)
        threshold = threshold_scaled * nfix / nfree
        return threshold

    def plot_ci(self, para):
        assert para in self.p_names, 'para must be one of ' + str(self.p_names)
        f, ax = plt.subplots()
        xx = self.trace_dict[para]['value']
        yy = self.trace_dict[para]['dchi']
        t = self.trace_dict[para]['threshold']
        spl = sp.interpolate.UnivariateSpline(xx, yy, k=2, s=0)
        allx = np.linspace(xx[0], xx[-1], 20000)
        ax.plot(xx, yy, '+')
        ax.plot(allx, spl(allx), '-', lw=1)
        ax.axhline(t, color='k', ls='--', lw=0.5)
        ax.axvline(self.params[para].value, color='k', ls='-', lw=0.5)
        ax.axvspan(*self.ci_values[para], alpha=0.1, color='b')
        if self.log:
            ax.semilogx()
        ax.set_xlabel('Parameter value')
        ax.set_ylabel(r'$\chi^2\left/\chi^2_0\right. - 1$')
        ax.set_title(para)


def conf_interval(
    minimizer,
    result,
    p_names=None,
    prob=0.95,
    limits=0.5,
    log=False,
    points=11,
    return_CIclass=False,
):
    """
    Calculate the confidence interval (CI) for parameters.

    The parameter for which the CI is calculated will be varied, while the
    remaining parameters are re-optimized to minimize the chi-square. The
    resulting chi-square is used to calculate the probability with a given
    statistic, i.e. chi-squared test.

    Parameters
    ----------
    minimizer : Minimizer
        The minimizer to use, holding objective function.
    result : MinimizerResult
        The result of running minimize().
    p_names : list, optional
        Names of the parameters for which the CI is calculated. If None
        (default), the CI is calculated for every parameter.
    prob : float, optional
        The probability for the confidence interval (<1). If None,
        the default is 0.95 (95 % confidence interval).
    limits : float, optional
        The limits (as a fraction of the original parameter value) within which
        to vary the parameters for identifiability analysis (default is 0.5).
        If ``log=False``, the parameter is varied from (1 - limits)*p to
        (1 + limits)*p, where p is the original value.
        If ``log=True``, the parameter is varied from p*limits to p/limits.
    log : bool, optional
        Whether to vary the parameter in a log (True) or a linear (False,
        default) scale.
    points : int, optional
        The number of points for which to calculate the profile likelihood over
        the given parameter range.
    return_CIclass : bool, optional
        When true, return the instantiated ``ConfidenceInterval`` class to
        access its methods directly (default=False).

    Returns
    -------
    output : dict
        A dictionary containing a list of ``(lower, upper)``-tuples containing
        the confidence bounds for each parameter.
    ci : ``ConfidenceInterval`` instance, optional
        Instantiated ``ConfidenceInterval`` class to access the attached methods.
    """
    assert (limits > 0) & (limits < 1), 'Please select a limits value between 0 and 1.'
    ci = ConfidenceInterval(minimizer, result, p_names, prob, log)
    output = ci.calc_all_ci(limits, points)
    if return_CIclass:
        return output, ci
    return output
