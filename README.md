# Proximal_Conditional_Gradient_Penalty_Method

This is a Matlab package that implements a single-loop proximal-conditional-gradient penalty method for solving compressed sensing (CS) instances with generalized Gaussian measurement noise:

$$
\min\limits_{x} ~ \|x\|_1\quad \text{s.t.}~\|y\|_p\le\sigma, ~\|x\|_\infty\le\|A^\dagger b\|_1+1, ~y=Ax-b,
$$


In every iteration of our algorithm, we apply one step of proximal gradient algorithm and one step of the conditional gradient algorithm for the following penalty function of the problem:

$$
\widetilde F_{\beta_k} = \|x\|_1+\delta_{\|\cdot\|\le \|A^\dagger b\|_1+1}+\delta_{\|\cdot\|_p\le \sigma}+\frac{\beta_k}{2}\|Ax-y-b\|^2,
$$

where $k\ge 0$ is the iteration number and $ \beta_{k}= \beta_0\sqrt{k+1} $ with $\beta_0>0$.
<br />


### Runcodes of the numerical experiments:

**plotBeta**: test the algorithm with $\beta_0\in \{0.1, 1, 10, 20, 50\}$ for randomly generated instances; plot the relative primal-dual gaps (see the definition in the paper) and the quantity $(\|Ax_{\rm out}-b\|_p-\sigma)_+$ at $x_{\rm out}$ returned by the algorithm <br />
**plotRecovery**: test the algorithm for a fixed random instance and $\beta_0=20$; plot the sequence $\{\|Ax^k-y^k-b\|\}$ and the comparison between the original signal and $x_{\rm out}$ returned by the algorithm <br />


### Other Matlab source codes are:


**proxFW**: Main function implementing of our proximal-condition-gradient penalty method <br />
**genGauss**: Subroutine for simulating random variables from the generalized Gaussian distribution . The implementation refers to [https://blogs.sas.com/content/iml/2016/09/21/simulate-generalized-gaussian-sas.html](https://blogs.sas.com/content/iml/2016/09/21/simulate-generalized-gaussian-sas.html)<br />
**aboxplot, quartile, colorgrad**: Subroutines for boxplots in the runcode **plotBeta** <br />. These codes are from [https://github.com/lemonzi/matlab/tree/master/aboxplot](https://github.com/lemonzi/matlab/tree/master/aboxplot)
**savefig**: Subroutine for saving figures
