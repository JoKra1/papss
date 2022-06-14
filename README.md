# papss: Penalized Additive Pupil Spline Solver

## Description:
The code in this repository allows to estimate penalized versions of the additive pupil model originally proposed by Hoeks & Levelt (1993) and then refined by Wierda et al. (2012) and later Denison et al. (2020) (for a review about applications see Fink et al. (2021)). More specifically, the model relies on the pupil response function recovered by Hoeks & Levelt (1993) as a basis function for an additive model (e.g., Wood, 2017) of pupil dilation over time. Each pupil basis is associated with a weight term, reflecting the magnitude of a spike in cognitive demand (Hoeks & Levelt, 1993; Wierda et al., 2012). ``papss`` also allows to add slope and/or intercept terms to the pupil dilation model (based on the refinements by Wierda et al. 2012 and Denison et al. 2020). As is typically done in **penalized additive models**, a penalty is imposed on these terms (see Eilers & Marx, 2010 or Wood, 2017 for a discussion of the penalties used here and Eilers & Marx, 1996 for a discussion of how penalties imposed directly on the coefficients should be interpreted). This enables the model to better recover continuous changes in demand compared to what would be possible without the penalty.

Demand weights are recovered by obtaining the non-negative least squares estimate (NNLS, see Bolduc et al., 2017) that minimizes the difference between the weighted sum of the pupil-response basis functions (i.e., our pupil-spline) and the observed pupil dilation time-course. NNLS estimation is performed because Hoeks & Levelt (1993) argued that there is no easy physical justification for negative spikes (i.e., these would imply that a decrease in demand, relative to baseline demand, leads to *active constriction* of the pupil). Solving the NNLS problem is achieved here by relying on projected gradient descent (PGD, see Ang, 2020a; Ang, 2020b; Bolduc, 2017), where after each (accelerated, see Sutskever et al. 2013) gradient step, the spike weights are constrained to be >= 0.

To estimate the optimal degree of spike penalization we rely on the generalized Fellner-Schall update for the smoothness penalty described by Wood & Fasiolo (2017). The update rule allows to maximize the restricted likelihood of the additive pupil model (i.e., pupil dilation = pupil spline + N(0, sigma^2)), after the optimal estimates of the spike weights have been obtained for a given degree of penalization (see Wood & Fasiolo, 2017 for details). The update equation was designed for unconstrained least squares problems, and is based on the default Hessian matrix for these models (see Wood, 2011). In practice, using the latter for this problem works reasonably well. Alternatively, ``papss`` offers to approximate the Hessian from projected gradient information using the BFGS update (Fletcher, 2000; Kim et al., 2006).

## Implementation details:
The definition of the model including model matrix setup is handled in R (R core team, 2021). Solving the penalized NNLS problem and optimizing the additive model is implemented in Eigen (Guennebaud et al., 2010) in C++. RCPP and RCPPEIGEN  (Bates & Eddelbuettel) act as interface between R and C++. 

## Installation:

### Windows-specific:
On Windows you first need to install RTools to be able to build this package. Instructions can be found [here](https://cran.r-project.org/bin/windows/Rtools/rtools40.html). Make sure that you add RTools to the path, as described in the section "Putting RTools" on the Path. As recommended in the article, you can use the following command to make sure everything was installed correctly:

```
Sys.which("make")
```

You can also use this command in case you are unsure whether you have already installed RTools in the past.

### macOS-specific
On macOS you will first need to install [Xcode](https://apps.apple.com/de/app/xcode/id497799835?mt=12) to get compiler support for R.

### OS-independent
Subsequently, you should be able to install the package right from this repository using functionality offered by the ``remotes`` package, as described for example [here](https://cran.r-project.org/web/packages/remotes/readme/README.html). Importantly, if you want to build the example vignette, you should set the ``build_vignettes`` argument to *True*. In that case you also need to set the ``dependencies`` argument to *True*, to install the optional dependencies needed for the vignettes (e.g. ``knitr`` and ``rmarkdown``):

```
install.packages("remotes")
remotes::install_github("JoKra1/papss",build_vignettes = T, dependencies = T)
```

Once the installation has been completed, you can inspect the vignettes using one of the following commands:

```
browseVignettes() # Search for 'papss'
vignette("artificial_data_analysis", package="papss") # Browse introduction vignette directly
```

There are two vignettes included at the moment. The *artificial_data_analysis* vignette is an introduction to how the
package can be used in practice. The *convergence_analysis* vignette shows how reliable the algorithm recovers a given solution from
different starting points.

## ToDo:
- Add CV score optimization for lambda?
- Implement difference penalty (Eilers & Marx, 1996; Wood, 2017)

## References:

1. Ang, A. (2020a). Accelerated gradient descent for large-scale optimization: On Gradient descent solving Quadratic problems—Guest lecture of MARO 201—Advanced Optimization. https://angms.science/doc/teaching/GDLS.pdf

2. Ang, A. (2020b). Nonnegative Least Squares—PGD, accelerated PGD and with restarts. https://angms.science/doc/NMF/nnls_pgd.pdf

3. Bates, D., & Eddelbuettel, D. (2013). Fast and Elegant Numerical Linear Algebra Using the RcppEigen Package. Journal of Statistical Software, 52(5). https://doi.org/10.18637/jss.v052.i05

4. Bolduc, E., Knee, G. C., Gauger, E. M., & Leach, J. (2017). Projected gradient descent algorithms for quantum state tomography. Npj Quantum Information, 3(1), 1–9. https://doi.org/10.1038/s41534-017-0043-1

5. Denison, R. N., Parker, J. A., & Carrasco, M. (2020). Modeling pupil responses to rapid sequential events. Behavior Research Methods, 52(5), 1991–2007. https://doi.org/10.3758/s13428-020-01368-6

6. Eilers, P. H. C., & Marx, B. D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science, 11(2), 89–121. https://doi.org/10.1214/ss/1038425655

7. Eilers, P., & Marx, B. (2010). Splines, knots, and penalties. https://doi.org/10.1002/WICS.125

8. Fink, L., Simola, J., Tavano, A., Lange, E. B., Wallot, S., & Laeng, B. (2021). From pre-processing to advanced dynamic modeling of pupil data. PsyArXiv. https://doi.org/10.31234/osf.io/wqvue

9. Fletcher, R. (2000). Practical Methods of Optimization. John Wiley & Sons, Incorporated. http://ebookcentral.proquest.com/lib/rug/detail.action?docID=1212544

10. Guennebaud, G., Jacob, B., & others. (2010). Eigen v3. http://eigen.tuxfamily.org

11. Hoeks, B., & Levelt, W. (1993). Pupillary dilation as a measure of attention: A quantitative system analysis. Behav. Res. Meth. Ins. C., 25, 16–26.

12. Kim, D., Sra, S., & Dhillon, I. S. (2006). A New Projected Quasi-Newton Approach for the Non-negative Least Squares Problem.

13. R Core Team. (2021). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing. https://www.R-project.org/

14. Sutskever, I., Martens, J., Dahl, G., & Hinton, G. (2013). On the importance of initialization and momentum in deep learning. Proceedings of the 30th International Conference on Machine Learning, 1139–1147. https://proceedings.mlr.press/v28/sutskever13.html

15. Tseng, P. (2008). On Accelerated Proximal Gradient Methods for Convex-Concave Optimization. https://www.mit.edu/~dimitrib/PTseng/papers.html

16. Wierda, S. M., van Rijn, H., Taatgen, N. A., & Martens, S. (2012). Pupil dilation deconvolution reveals the dynamics of attention at high temporal resolution. Proceedings of the National Academy of Sciences of the United States of America, 109(22), 8456–8460. https://doi.org/10.1073/pnas.1201858109

17. Wood, S. N., & Fasiolo, M. (2017). A generalized Fellner-Schall method for smoothing parameter optimization with application to Tweedie location, scale and shape models. Biometrics, 73(4), 1071–1081. https://doi.org/10.1111/biom.12666

18. Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models: Estimation of Semiparametric Generalized Linear Models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(1), 3–36. https://doi.org/10.1111/j.1467-9868.2010.00749.x

19. Wood, S. N. (2017). Generalized Additive Models: An Introduction with R, Second Edition (2nd ed.). Chapman and Hall/CRC.
