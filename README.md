# identifiability - Parameter identifiability analysis in Python

This module performs parameter identifiability
analysis to calculate and plot confidence intervals based on a profile-likelihood. 
The code is adapted from [LMFIT](https://lmfit.github.io/lmfit-py/), with custom
functions to select the range for parameter scanning and for plotting the profile 
likelihood. The significance is assessed with the chi-squared distribution. 
Optimization runs can be performed in parallel (using the `multiprocessing` module).

## Installation

`identifiability` is a pure-Python module. The latest development version can be 
installed with
```bash
$ pip install https://github.com/jmrohwer/identifiability/archive/refs/heads/main.zip
```

The latest stable release is available on PyPI:
```bash
$ pip install identifiability
```
The module can be used in combination with [PySCeS](https://pysces.github.io) for 
simulation and parameter estimation of kinetic models using the `CVODE` solver. When 
performing identifiability analysis in parallel using `multiprocessing`, additional 
dependences are required; these can be installed with:
```bash
$ pip install "identifiability[pyscesmp]"
```

## API change

> **NOTE:**  
>The API of the `identifiability` module has changed since version 0.5.0. The `conf_interval()` 
helper function has been removed due to incompatibilities with Python 3.13. The `ConfidenceInterval`
class now has to be instantiated directly.

## Basic usage

For background, the reader is referred to the section on *Calculation of confidence 
intervals* in the [LMFIT documentation](https://lmfit.github.io/lmfit-py/confidence.html).

To start the identifiability analysis, the user first needs to have performed a 
parameter estimation with LMFIT. The method for estimating confidence intervals 
takes an instantiated LMFIT 
[Minimizer](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.Minimizer)
object and a 
[MinimizerResult](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.MinimizerResult)
object as input.

A typical workflow would entail:
```python
>>> from identifiability import ConfidenceInterval
>>> ci = ConfidenceInterval(mini, result)
>>> ci.calc_all_ci()    # returns OrderedDict of parameter names and 
                        # corresponding confidence intervals 
>>> ci.plot_ci('a')     # plots confidence interval for parameter 'a'
>>> ci.plot_all_ci()    # plots confidence intervals for all parameters
```

When using the [Model](https://lmfit.github.io/lmfit-py/model.html) wrapper of LMFIT 
to perform the parameter estimation and model fit, the instantiated 
[ModelResult](https://lmfit.github.io/lmfit-py/model.html#lmfit.model.ModelResult)
object should be passed twice to the `conf_interval()` method, instead of the
[Minimizer](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.Minimizer)
and 
[MinimizerResult](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.MinimizerResult)
(see above). In this case the function call would be:
```python
>>> ci = ConfidenceInterval(modelresult, modelresult)
```

Once a profile likelihood has been calculated, the same data can be used to calculate 
the confidence interval for a different probability, thus avoiding the 
computationally intensive re-calculation of the profile likelihood:

```python
>>> ci.calc_all_ci(prob=0.8)
```

The method for calculating the confidence intervals (`calc_all_ci()`) has several
additional options to specify the parameters to be varied, the probability for the
confidence interval, the limits for the parameter variation, and whether to use a
linear or logarithmic scale for the parameter variation. See below.

### Docstrings of the `ConfidenceInterval` class and selected methods

```python
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

    def plot_all_ci(self):
        """Plot profile likelihoods with confidence intervals for all parameters."""
        
    def plot_ci(self, para, ax=None):
        """Plot profile likelihood with confidence interval for single parameter.

        Parameters
        ----------
        para : str
            The parameter name for which to plot the profile likelihood.
        ax: matplotlib.axes.Axes, optional
            Matplotlib Axes object to plot on. If None, a new figure and axes will be created.
        """


```


© Johann M. Rohwer, 2023&ndash;2025