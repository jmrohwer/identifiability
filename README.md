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

## Basic usage

For background, the reader is referred to the section on *Calculation of confidence 
intervals* in the [LMFIT documentation](https://lmfit.github.io/lmfit-py/confidence.html).

To start the indentifiability analysis, the user first needs to have performed a 
parameter estimation with LMFIT. The method for estimating confidence intervals 
takes an instantiated LMFIT 
[Minimizer](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.Minimizer)
object and a 
[MinimizerResult](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.MinimizerResult)
object as input.

A typical workflow would entail:
```python
>>> from identifiability import conf_interval
>>> c = conf_interval(
        mini, result, prob=0.95, limits=0.5, log=False, points=11, return_CIclass=True
    )
>>> print(c[0])  # OrderedDict of parameter names and corresponding confidence intervals
>>> c[1].plot_ci('a')   # plot confidence interval for parameter 'a'
>>> c[1].plot_all_ci()  # plot confidence intervals for all parameters
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
>>> c = conf_interval(
        modelresult, modelresult, prob=0.95, limits=0.5, 
        log=False, points=11, return_CIclass=True
    )
```

Once a profile likelihood has been calculated, the same data can be used to calculate 
the confidence interval for a different probability, thus avoiding the 
computationally intensive re-calculation of the profile likelihood:

```python
>>> c[1].calc_all_ci(prob=0.8)
```

### Docstring of the `conf_interval` method

```python
def conf_interval(
    minimizer,
    result,
    p_names=None,
    prob=0.95,
    limits=0.5,
    log=False,
    points=11,
    method='leastsq',
    return_CIclass=False,
    mp=True,
):
    """
    Calculate the confidence interval (CI) for parameters.

    The parameter for which the CI is calculated will be varied, while the
    remaining parameters are re-optimized to minimize the chi-square. The
    resulting chi-square is used to calculate the probability with a given
    statistic, i.e. chi-squared test.

    Parameters
    ----------
    minimizer : Minimizer or ModelResult
        The minimizer to use, holding objective function.
    result : MinimizerResult or ModelResult
        The result of running Minimizer.minimize() or Model.fit().
    p_names : list, optional
        Names of the parameters for which the CI is calculated. If None
        (default), the CI is calculated for every parameter.
    prob : float, optional
        The probability for the confidence interval (<1). If None,
        the default is 0.95 (95 % confidence interval).
    limits : float, optional
        The limits (as a fraction of the original parameter value) within which
        to vary the parameters for identifiability analysis (default is 0.5).
        If ``log=False``, the parameter is varied from p*limits to p*(2 - limits), 
        where p is the original value.
        If ``log=True``, the parameter is varied from p*limits to p/limits.
    log : bool, optional
        Whether to vary the parameter in a log (True) or a linear (False,
        default) scale.
    points : int, optional
        The number of points for which to calculate the profile likelihood over
        the given parameter range.
    method : str, optional
        The lmfit mimimize() method to use (default='leastsq')
    return_CIclass : bool, optional
        When true, return the instantiated ``ConfidenceInterval`` class to
        access its methods directly (default=False).
    mp : bool, optional
        Run the optimization in parallel using ``multiprocessing`` (default=True)

    Returns
    -------
    output : dict
        A dictionary containing a list of ``(lower, upper)``-tuples containing
        the confidence bounds for each parameter.
    ci : ``ConfidenceInterval`` instance, optional
        Instantiated ``ConfidenceInterval`` class to access the attached methods.
    """
```

© Johann M. Rohwer, 2022