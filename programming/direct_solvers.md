# LU factorization
Let A be a general, square, invertible matrix: 

<a href="https://www.codecogs.com/eqnedit.php?latex=A&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;\end{bmatrix}" title="A = \begin{bmatrix} 9 & 0 & 3 & 0 \\ 0 & 8 & 0 & 0 \\ 0 & 2 & 6 & 0 \\ 1 & 0 & 0 & 5 \end{bmatrix}" /></a>, 

Calculate its LUP decomposition `A = P*L*U`: 
```csharp
Matrix A = Matrix.CreateFromArray(new double[,] 
{
    { 9.0, 0.0, 3.0, 0.0 },
    { 0.0, 8.0, 0.0, 0.0 },
    { 0.0, 2.0, 6.0, 0.0 },
    { 1.0, 0.0, 0.0, 5.0 }
});

LUFactorization lu = A.FactorLU(inPlace: true);
//Vector diag = A.GetDiagonal(); // A is now overwritten. This will throw an exception.
```

The parameter `inPlace` specifies whether the `LUFactorization` object should use the same array as the original `Matrix` object and overwrite it with the factors L, U. This would avoid having to allocate new memory and could be necessary if the matrix is too large to be stored twice. However, by doing so the original `Matrix` object should no longer be used for matrix operations. This behavior will be present in most factorizations henceforth.

## Solving linear systems
Once the LUP factorization is computed, it can be used to solve one or more linear systems:
```csharp
LUFactorization lu;
Vector b = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
Vector x = Vector.CreateZero(4);
lu.SolveLinearSystem(b, x);
```

## Determinant and inverse
The LUP factorization can also be used to calculate the determinant and the inverse matrix, in case the latter is explicitly needed:
```csharp
LUFactorization lu;
double det = lu.CalcDeterminant();
Matrix invA = lu.Invert(inPlace: true);
```

## Sparse matrix formats
LU(P) decomposition is also supported for CSC storage format. The fill-in due to the factorization cannot be contained within the original CSC matrix and new memory always needs to be allocated.

```csharp
double[] valuesA = { 9, 1, 8, 2, 3, 6, 5 };
int[] rowIndicesA = { 0, 3, 1, 2, 0 };
int[] colOffsetsA = { 0, 2, 4, 6, 7 };
CscMatrix A = CscMatrix.CreateFromArrays(4, 4, valuesA, rowIndicesA, colOffsetsA, checkInput: true);

// Factorize and solve linear system
LUCSparseNet lu = LUCSparseNet.Factorize(A);
Vector b = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
Vector x = Vector.CreateZero(4);
lu.SolveLinearSystem(b, x);
```

# Cholesky / LDL factorization
For a symmetric positive definite matrix 

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\begin{bmatrix}&space;10&space;&&space;0&space;&&space;2&space;&&space;0\\&space;0&space;&&space;10&space;&&space;0&space;&&space;2\\&space;2&space;&&space;0&space;&&space;10&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;0&space;&&space;10&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\begin{bmatrix}&space;10&space;&&space;0&space;&&space;2&space;&&space;0\\&space;0&space;&&space;10&space;&&space;0&space;&&space;2\\&space;2&space;&&space;0&space;&&space;10&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;0&space;&&space;10&space;\end{bmatrix}" title="A=\begin{bmatrix} 10 & 0 & 2 & 0\\ 0 & 10 & 0 & 2\\ 2 & 0 & 10 & 0 \\ 0 & 2 & 0 & 10 \end{bmatrix}" /></a>

we can calculate its Cholesky or LDL decomposition instead of LUP. In this case the packed matrix format of `SymmetricMatrix` can also be used to save up memory:
```csharp
Matrix A = Matrix.CreateFromArray(new double[,] 
{
    { 10.0, 0.0, 2.0, 0.0 },
    { 0.0, 10.0, 0.0, 2.0 },
    { 2.0, 0.0, 10.0, 0.0 },
    { 0.0, 2.0, 0.0, 10.0 }
});
SymmetricMatrix symA = SymmetricMatrix.CreateFromMatrix(A);

// Factorize
CholeskyFull ll = A.FactorCholesky(inPlace: true);
CholeskyPacked symLL = symA.FactorCholesky();

// Solve linear system
Vector b = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
Vector x = Vector.CreateZero(4);
ll.SolveLinearSystem(b, x);
symLL.SolveLinearSystem(b, x);

// Determinant, inverse 
double det = ll.CalcDeterminant();
det = symLL.CalcDeterminant();
Matrix invA = ll.Invert(inPlace: true);
SymmetricMatrix symInvA = symLL.Invert(inPlace: true);
```

## Sparse matrix formats
Cholesky/LDL factorizations are also supported for Skyline and CSC sparse matrix formats. In the case of CSC format, only the upper triangle is stored and new memory is always allocated for the factorization. On the other hand, the Skyline format stores exactly the entries that will be filled during factorization and can thus be overwritten by the factorized data.
```csharp
// Skyline format: it contains the fill-in due to factorization
double[] skyValuesA = { 10.0, 10.0, 10.0, 0.0, 2.0, 10.0, 0.0, 2.0 };
int[] skyDiagOffsetsA = { 0, 1, 2, 4, 6 };
SkylineMatrix skyA = SkylineMatrix.CreateFromArrays(4, skyValuesA, skyDiagOffsetsA, checkInput: true);

// CSC format: only upper triangle
double[] cscValuesA = { 10.0, 10.0, 2.0, 10.0, 2.0, 10.0 };
int[] cscRowIndicesA = { 0, 1, 0, 2, 1, 3 };
int[] cscColOffsetsA = { 0, 1, 2, 4, 6 };
SymmetricCscMatrix cscA = SymmetricCscMatrix.CreateFromArrays(4, cscValuesA, cscRowIndicesA, cscColOffsetsA, checkInput: true);

// Factorize
LdlSkyline skyLDL = skyA.FactorLdl(inPlace: true); // Overwrites the entries of A with the factorization.
CholeskyCSparseNet cscLL = CholeskyCSparseNet.Factorize(cscA); // Allocates new memory for the factorization.

// Solve linear system
Vector b = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
Vector x = Vector.CreateZero(4);
cscLL.SolveLinearSystem(b, x);
skyLDL.SolveLinearSystem(b, x);

// Solve multiple linear systems at once
Matrix B = Matrix.CreateFromArray(new double[,] 
{
    { 1.0, 5.0 },
    { 2.0, 6.0 },
    { 3.0, 7.0 },
    { 4.0, 8.0 }
}); 
Matrix X = Matrix.CreateZero(4, 2);
skyLDL.SolveLinearSystems(B, X);
```

# QR and LQ factorization
Let A, B be two general matrices (rectangular or square):

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\begin{bmatrix}&space;1&space;&&space;5&space;&&space;9\\&space;2&space;&&space;6&space;&&space;10\\&space;3&space;&&space;7&space;&&space;11\\&space;4&space;&&space;8&space;&&space;12&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\begin{bmatrix}&space;1&space;&&space;5&space;&&space;9\\&space;2&space;&&space;6&space;&&space;10\\&space;3&space;&&space;7&space;&&space;11\\&space;4&space;&&space;8&space;&&space;12&space;\end{bmatrix}" title="A=\begin{bmatrix} 1 & 5 & 9\\ 2 & 6 & 10\\ 3 & 7 & 11\\ 4 & 8 & 12 \end{bmatrix}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=B=\begin{bmatrix}&space;1&space;&&space;4&space;&&space;7&space;&&space;10\\&space;2&space;&&space;5&space;&&space;8&space;&&space;11\\&space;3&space;&&space;6&space;&&space;9&space;&&space;12&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?B=\begin{bmatrix}&space;1&space;&&space;4&space;&&space;7&space;&&space;10\\&space;2&space;&&space;5&space;&&space;8&space;&&space;11\\&space;3&space;&&space;6&space;&&space;9&space;&&space;12&space;\end{bmatrix}" title="B=\begin{bmatrix} 1 & 4 & 7 & 10\\ 2 & 5 & 8 & 11\\ 3 & 6 & 9 & 12 \end{bmatrix}" /></a>

We can calculate the QR and LQ decompositions for A and B respectively. Then these can be used to solve a least-squares and minimum-norm problem respectively:

```csharp
Matrix A = Matrix.CreateFromArray(new double[,] 
{
    { 1.0, 5.0,  9.0 },
    { 2.0, 6.0, 10.0 },
    { 3.0, 7.0, 11.0 },
    { 4.0, 8.0, 12.0 }
});
Matrix B = Matrix.CreateFromArray(new double[,] 
{
    { 1.0, 4.0, 7.0, 10.0 },
    { 2.0, 5.0, 8.0, 11.0 },
    { 3.0, 6.0, 9.0, 12.0 },
});

// Factorize
bool inPlace = true;
QRFactorization qr = A.FactorQR(inPlace);
LQFactorization lq = A.FactorLQ(inPlace);

// Solve least-squares problem
Vector b1 = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
Vector x1 = qr.SolveLeastSquares(b1);

// Solve min-norm problem
Vector b2 = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0 });
Vector x2 = lq.SolveMinNorm(b2);
```

# Reordering
When applying LU, Cholesky and LDL factorizations for sparse matrices, fill-in occurs: a lot of zero entries are filled with non-zero values. This increases the memory and computational time requirements. The amount of fill-in depends on the matrix bandwidth and can be reduced by reordering the rows and columns of the matrix before factorization. The following example demonstrates a simplification of the steps taken to use a direct sparse solver for the linear system resulting by the Finite Element Method or a similar one: 

In FEM the rows and columns of the linear system's matrix (called stiffness matrix) represent freedom degrees at nodes of the discretized domain. Having a non-zero entry A<sub>ij</sub> means that there is coupling between the freedom degrees i, j. This coupling is the result of a finite element connececting the nodes where these freedom degrees are defined. When we reorder the rows and columns of the stiffness matrix, we also change the order used to number the freedom degrees. Therefore the reordering will be also used to update data of the finite element model. Assume that there are 3 finite elements and the stiffness matrix of each one is equal to:

<a href="https://www.codecogs.com/eqnedit.php?latex=k&space;=&space;\begin{bmatrix}&space;2&space;&&space;1\\&space;1&space;&&space;2&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k&space;=&space;\begin{bmatrix}&space;2&space;&&space;1\\&space;1&space;&&space;2&space;\end{bmatrix}" title="k = \begin{bmatrix} 2 & 1\\ 1 & 2 \end{bmatrix}" /></a>

Also assume that due to the original order of freedom degrees, the global (for the whole domain) stiffness matrix will be:

<a href="https://www.codecogs.com/eqnedit.php?latex=K&space;=&space;\begin{bmatrix}&space;2&space;&&space;&&space;1&space;&&space;\\&space;&&space;4&space;&&space;1&space;&&space;1\\&space;1&space;&&space;1&space;&&space;4&space;&&space;\\&space;&&space;1&space;&&space;&&space;2&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K&space;=&space;\begin{bmatrix}&space;2&space;&&space;&&space;1&space;&&space;\\&space;&&space;4&space;&&space;1&space;&&space;1\\&space;1&space;&&space;1&space;&&space;4&space;&&space;\\&space;&&space;1&space;&&space;&&space;2&space;\end{bmatrix}" title="K = \begin{bmatrix} 2 & & 1 & \\ & 4 & 1 & 1\\ 1 & 1 & 4 & \\ & 1 & & 2 \end{bmatrix}" /></a>

We want to find a permutation of the rows and columns that minimizes the bandwidth and thus the fill-in of the global stiffness matrix K. To do so, we only need to know the connectivity between the finite elements and their nodes. From this data we can construct arrays that to which K<sub>ij</sub> we must add each k<sub>ij</sub>. For the matrix K shown above these would be:

```csharp
// Mapping arrays. In this case rows and columns are identical, since all matrices is symmetric.
var localIndices = new List<int[]>();
for (int i = 0; i < 3; ++i) localIndices[i] = new int[] { 0, 1 };
var globalIndices = new List<int[]>();
globalIndices[0] = new int[] { 0, 2 };
globalIndices[1] = new int[] { 2, 1 };
globalIndices[2] = new int[] { 1, 3 };
```
We can now determine the sparsity pattern of the global stiffness matrix K, without building it. Note that the element stiffness matrices k do not to be calculated at this point. The mapping arrays above are determined by the element connectivity. They will only be needed when assembling the final global stiffness matrix, after the reordering is finished.

```csharp
SparsityPatternSymmetric pattern = SparsityPatternSymmetric.CreateEmpty(4);
for (int i = 0; i < 3; ++i) pattern.ConnectIndices(globalIndices[i]);
```

 After determining the sparsity pattern we can find a fill-reducing permutation. Here we will use the Approximate Minimum Degree algorithm (AMD). This permutation can be used to reorder the rows and columns of the global stiffness matrix. To do this we essentially need to modify the global order of freedom degrees:

```csharp
(int[] permutation, bool oldToNew) = FindPermutation(SparsityPatternSymmetric pattern);

// Update the order of global freedom degrees. The following assumes that oldToNew = false
for (int i = 0; i < 3; ++i)
{
    int[] dofs = globalIndices[i];
    for (int j = 0; j < dofs.Length; ++j) dofs[j] = dofs[permutation[j]];
}
```

Now we can calculate the element stiffness matrices and assemble them into the reordered global stiffness matrix K. Here we will use the CSC format for the upper triangle of K:
```csharp
// Calculate element matrices. This could require a lot of work in FEM-like methods. Thankfully we only need to do this once, just before assembling the reordered global matrix.
var elementMatrices = new List<Matrix>();
for (int i = 0; i < 3; ++i)
{
    Matrix k = Matrix.CreateZero(2, 2);
    k[0, 0] = 2.0;
    k[0, 1] = 1.0;
    k[1, 0] = 1.0;
    k[1, 1] = 2.0;
    elementMatrices[i] = k;
}

// Assemble element matrices into the reordered global matrix
DokSymmetric dok = DokSymmetric.CreateEmpty(4);
for (int i = 0; i < 3; ++i)
{
    dokForSymCsc.AddSubmatrixSymmetric(elementMatrices[i], localIndices[i], globalIndices[i]);
}
```