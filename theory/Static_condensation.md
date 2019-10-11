# Schur complement
Let us seprarate the rows and columns of a symmetric matrix A into 2 groups and define the block submatrices: 

![condensation_decomposition_matrix](img/condensation_decomposition_matrix.png)

If A<sub>22</sub> is invertible, the Schur complement of the block A<sub>22</sub> is defined as:

![condensation_matrix](img/condensation_matrix.png)

The Schur complement is important in many applications, such as domain decomposition methods, homogenization, etc, and has some interesting properties: e.g. if A is symmetric positive definite, then so is the Schur complement of A<sub>22</sub> and vice-versa. For more details see this [wikipedia article](https://en.wikipedia.org/wiki/Schur_complement).


# Static condensation
Let us apply the previous decomposition to a linear system `A*x = b`:

![condensation_decomposition_matrix](img/condensation_decomposition_matrix.png) ![condensation_decomposition_vectors](img/condensation_decomposition_vectors.png)

Suppose we need to work only with the group 1 of rows/columns, without completely disregarding the contribution of group 2. In this case we typically condense the submatrix A<sub>22</sub> using its Schur complement. We also condense the subvector b<sub>2</sub>:

![condensation_rhs](img/condensation_rhs.png)

Then we can solve the system with respect to x<sub>1</sub> only:

![condensation_system](img/condensation_system.png)

If x<sub>2</sub> is needed, it can be calculated afterwards by solving the linear system:

![condensation_lhs](img/condensation_lhs.png)