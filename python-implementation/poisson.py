#from https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.stats.poisson.html#scipy.stats.poisson
from scipy.stats import poisson
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)

mu = 0.6

mean, var, skew, kurt = poisson.stats(mu, moments='mvsk')

x = np.arange(poisson.ppf(0.01, mu), poisson.ppf(0.99, mu))
ax.plot(x, poisson.pmf(x, mu), 'bo', ms=8, label='poisson pmf')
ax.vlines(x, 0, poisson.pmf(x, mu), colors='b', lw=5, alpha=0.5)

rv = poisson(mu)
ax.vlines(x, 0, rv.pmf(x), colors='k', linestyles='-', lw=1,label='frozen pmf')
ax.legend(loc='best', frameon=False)
plt.show()
