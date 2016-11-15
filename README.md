# Brief description

NumDiff provides a modern Fortran interface for computing the Jacobian (derivative) matrix of `m` nonlinear functions which depend on `n` variables. The Jacobian matrix is required for various applications, including numerical optimization. The library also provides for computing the sparsity of this matrix, and returning the Jacobian in sparse or dense form.

# Status

This is currently an experimental work in progress and is not production ready. The goal is a comprehensive library that contains a full suite of computationally efficient implementations of algorithms for sparsity determination and numerical differentiation.

#### To Do:

- [ ] General code cleanup
- [x] Computing the nonlinear sparsity pattern
  - [x] Specified by the user
  - [x] Assume all elements `true`
  - [x] Three random points within variable bounds
- [ ] Computing the linear sparsity pattern
- [ ] Add additional finite different gradient methods
  - [x] Forward differences
  - [x] Central differences
  - [x] Backward differences
  - [ ] Higher-order differences...
- [x] Perturbations should respect variable bounds
- [x] Ability to use different methods for different columns
- [ ] Jacobian partitioning to compute multiple columns at the same time
- [ ] Add other gradient methods?
- [ ] Also compute Hessian matrix?
- [ ] OpenMP or Coarrays for parallelization
- [ ] Testing for computational efficiency

# License

The NumDiff source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/NumDiff/blob/master/LICENSE) (BSD-style).
