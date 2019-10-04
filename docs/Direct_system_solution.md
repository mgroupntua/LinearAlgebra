# LU factorization
To solve a linear system `A*x = b`, the square matrix `A (n x n)` must first be decomposed into simpler matrices that we can work with easier. If we record the operations done during [Gauss elimination](https://en.wikipedia.org/wiki/Gaussian_elimination) in matrices, then the LU factorization is obtained `A = L*U`, where `L (n x n)` is a lower triangular matrix with its diagonal entries being 1 and `U (n x n)` is an upper triangular matrix. To solve the linear system after factorizing the matrix, we perform [back substitution](https://algowiki-project.org/en/Backward_substitution) and [forward substitution](https://en.wikipedia.org/wiki/Triangular_matrix#Forward_and_back_substitution) with its factors L, U respectively: 
```
L * y = b (back substitution)
U * x = y (forward substitution)
```

Only matrices that do not need pivoting during Gauss elimination admit a LU factorization. In general pivoting is needed, which leads to the LUP factorization `A = P*L*U` where `P (n x n)` is a permutation matrix. All matrices admit an LUP factorization, but if the matrix is singular, the last diagonal entries of U will be 0 and the factorization cannot be used to solve linear systems. To solve the linear system, after the factorization, we perform again back and forward substitutions, but also multiply with the inverse of P, which is its transpose:
```
z = P^T * b (permutation)
L * y = z (back substitution)
U * x = y (forward substitution)
```

For more details about LU and LUP factorizations, see [Numerical Linear Algebra - Tefethen, Bau](https://books.google.gr/books/about/Numerical_Linear_Algebra.html?id=JaPtxOytY7kC&redir_esc=y) or [this wikipedia article](https://en.wikipedia.org/wiki/LU_decomposition).

# Cholesky factorization
If the linear system matrix A is (symmetric) positive definite, then it will always admit a Cholesky factorization `A=L*L^T`, where L is a lower triangular matrix. Cholesky factorization performs half the operations needed by LU, does not need pivoting and is always the better choice, provided the matrix is positive definite. After factorizing the matrix, we solve the linear system:
```
L * y = b (back substitution)
L^T * x = y (forward substitution)
```

A variant of Cholesky factorization is the LDL factorization  `A = L*D*L^T`, where D is a diagonal matrix and L is a lower triangular matrix with its diagonal entries being 1. The LDL factorization avoids calculating square roots and can be used for some symmetric indefinite matrices, for which Cholesky is unsuitable. To solve the linear system after factorizing the matrix, we must multiply with the inverse of D, which is trivial to calculate: D<sup>-1</sup><sub>ii</sub> = 1 / D<sub>ii</sub>.
```
L * z = b (back substitution)
y = inv(D) * z (diagonal inversion)
L^T * x = y (forward substitution)
```

For more details about Cholesky and LDL factorizations, see [Numerical Linear Algebra - Tefethen, Bau](https://books.google.gr/books/about/Numerical_Linear_Algebra.html?id=JaPtxOytY7kC&redir_esc=y) or [this wikipedia article](https://en.wikipedia.org/wiki/Cholesky_decomposition).

# Sparse matrices
If the linear sytem matrix A is sparse, then factorization will cause a lot of the zero entries to become non-zeros. This is called fill-in and it usually requires new data to be allocated to store the factorization matrices and extra computational work for processing the new non-zero entries. There are many iterative algorithms for solving linear systems that do not factorize the matrix and thus avoid fill-in. These are sometimes more efficient than the direct linear solution algorithms shown above or they can even be the only alternative if the sparse matrix is so large, that its factorization cannot be stored in the computer memory.

# Matrix inversion
Usually the inverse A<sup>-1</sup> of a matrix A is not calculated explicitly. Instead the matrix is factorized using LUP, Cholesky or LDL and then the resulting factors are used to solve one or more linear systems. However, if the inverse is indeed needed explicitly, then it can be obtained by the using the factorization data, e.g. by solving a set of linear systems L*U * X<sub>j</sub> = I<sub>j</sub>, where I<sub>j</sub> are the columns of the identity matrix and X<sub>j</sub> are the columns of A<sup>-1</sup>