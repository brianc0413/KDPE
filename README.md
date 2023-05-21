# KDPE
All code for the simulations presented in "Kernel Debiased Plug-in Estimation" is provided in this repository. To ensure that the code runs, our simulations require the packages $\texttt{SuperLearner}$, $\texttt{tmle}$, $\texttt{CVXR}$, and $\texttt{ggplot2}$. 

The main function to test our algorithm is $\texttt{kdpe}(\cdot)$, which is implemented in the helper file. In particular, it takes in the following arguments:

- $\texttt{Y}$: binary outcome vector
- $\texttt{A}$: binary treatment vector
- $\texttt{X}$: covariate / baseline characteristics vector
- $\texttt{g}$: a function $g(x)$ that serves as the initial estimate of $P(A=1|X=x)$, must work with predict function
- $\texttt{Q}$: a function $Q(a,x)$ that serves as the intial estimate of $P(Y=1|A=a, X=x)$, must work with predict function
- $\texttt{density_bound}$: $c$, the density bound for KDPE, as described in Section 4
- $\texttt{converge_tol}$: $\gamma$, the convergence tolerance, as described in Section 4

In this implementation, we use our convergence tolerance function to be $l(P,P')= \sum_{i=1}^n[P-P']^2$ (i.e. squared L-2 norm), and set $\texttt{density_bound} = 0.01, \ \texttt{converge_tol} = 0.001$ by default. 

The function returns a vector of length 4, which contains the following variables in this order:
- $\hat{\psi}_{\text{ATE}}$: ATE estimate
- $\hat{\psi}_{\text{RR}}$: RR estimate
- $\hat{\psi}_{\text{OR}}$: OR estimate
- $\texttt{iter_count}$: the number of iterations until convergence

This implementation of KDPE is meant to be a proof-of-concept, rather than a working version of this algorithm. In particular, it is not adapted for examples outside of the simulation set-up we provide in this paper. All code presented here is included for reproducibility of the results shown in "Kernel Debiased Plug-in Estimation."
