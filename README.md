# Brief description

NumDiff provides a modern Fortran interface for computing the Jacobian (derivative) matrix of `m` nonlinear functions which depend on `n` variables. The Jacobian matrix is required for various applications, including numerical optimization. The library also provides for computing the sparsity of this matrix, and returning the Jacobian in sparse or dense form.

# Status

This is currently an experimental work in progress and is not production ready. The goal is a comprehensive library that contains a full suite of computationally efficient implementations of algorithms for sparsity determination and numerical differentiation.

#### To Do:

- [x] Computing the nonlinear sparsity pattern
  - [x] Specified by the user
  - [x] Assume all elements `true`
  - [x] Three random points within variable bounds
- [x] Various order finite different gradient methods
  - [x] 2-point (backward, forward)
  - [x] 3-point (backward, central, forward)
  - [x] 4-point (backward 1, backward 2, forward 1, forward 2)
  - [x] 5-point (backward 1, backward 2, central, forward 1, forward 2)
  - [x] 6-point (backward 1, backward 2, backward 3, forward 1, forward 2, forward 3)
- [x] Perturbations should respect variable bounds
- [x] Neville's process
- [x] Ability to use different methods for different columns
- [x] Jacobian partitioning to compute multiple columns at the same time
- [ ] Estimate the optimal perturbation step size
- [ ] Computing the linear sparsity pattern (constant elements of Jacobian)
- [ ] Add other gradient methods?
- [ ] Also compute Hessian matrix?
- [ ] OpenMP or Coarrays for parallelization
- [ ] Testing for computational efficiency
- [ ] General code cleanup

# License

The NumDiff source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/NumDiff/blob/master/LICENSE) (BSD-style).
