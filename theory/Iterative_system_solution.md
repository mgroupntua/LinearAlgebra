# Introduction
Iterative methods can solve a linear system **A*x=b** by applying vector-vector operations and matrix-vector multiplications. This makes them very efficient for solving linear systems with large sparse matrices, which could not be factorized due to memory constraints. Furthermore they can be applied in cases where the linear system matrix is not explicitlty formed, such as matrix-free methods, subproblems of domain decomposition methods, etc.

# Conjugate Gradient (CG)
The Conjugate Gradient is the predominant iterative method for solving linear systems with a symmetric positive definite matrix. In exact arithmetic CG would converge in N iterations at  most, where N is the number of rows/columns of the matrix. In practice, this is not guaranteeed due to round-off errors and often we need faster convergence, therefore preconditioning is required. For more details see [Painless Conjugate Gradient](https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf).

# Minimum Residual (MINRES)
The Minimum Residual method is a variant of CG that is also applicable to symmetric indefinite linear systems. However it is not as efficient as CG. For more details see [Minimum Residual Method](http://mathworld.wolfram.com/MinimalResidualMethod.html).

# Generalized Minimum Residual (GMRES)
The Generalized Minimum Residual method is a generalization of MINRES, which be used to solve nonsymmetric linear systems as well. Its memory requirements and computational cost is higher than MINRES. For more details, see [Iterative Methods for Sparse Linear Systems](https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf).

# Preconditioning
Iterative methods converge to a solution after a number of iterations. This number depends on the condition number of the matrix, which can be defined as the ratio of maximum to minimum eigenvalue. In practice iterative methods are accompanied by a preconditioner, which transforms the linear system, such that the new matrix has a lower condition number and the iterative method converges faster. A preconditioner is an approximation of the inverse of the original linear system matrix. It should be cost efficient to calculate and apply and at the same time effective in reducing the number of iterations. There are general purpose and problem dependent preconditioners, e.g.:
* Jacobi or diagonal preconditioner. This preconditioner is obtained by inverting the diagonal of the linear system matrix. For more details see [Iterative Methods for Sparse Linear Systems](https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf).
* Dirichlet preconditioner for Finite Element Tearing and Interconnecting (FETI) domain decomposition methods. For more details see [FETI-DP: a dual-primal unified FETI method](https://web.stanford.edu/group/fpc/Publications/FETI.pdf).
