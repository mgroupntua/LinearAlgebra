In the next examples the right hand side vector will be:

 <a href="https://www.codecogs.com/eqnedit.php?latex=b=\begin{bmatrix}&space;1\\&space;2\\&space;3\\&space;4&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?b=\begin{bmatrix}&space;1\\&space;2\\&space;3\\&space;4&space;\end{bmatrix}" title="b=\begin{bmatrix} 1\\ 2\\ 3\\ 4 \end{bmatrix}" /></a>

# Conjugate Gradient
Let us solve the linear system `A * x = b`, where:

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\begin{bmatrix}&space;10&space;&&space;0&space;&&space;2&space;&&space;0\\&space;0&space;&&space;10&space;&&space;0&space;&&space;2\\&space;2&space;&&space;0&space;&&space;10&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;0&space;&&space;10&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\begin{bmatrix}&space;10&space;&&space;0&space;&&space;2&space;&&space;0\\&space;0&space;&&space;10&space;&&space;0&space;&&space;2\\&space;2&space;&&space;0&space;&&space;10&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;0&space;&&space;10&space;\end{bmatrix}" title="A=\begin{bmatrix} 10 & 0 & 2 & 0\\ 0 & 10 & 0 & 2\\ 2 & 0 & 10 & 0 \\ 0 & 2 & 0 & 10 \end{bmatrix}" /></a>   

A is symmetric positive definite. We will solve the linear system using the Conjugate Gradient method. The Jacobi preconditioner will be used, due to its simplicity. The matrix will be in CSR format, since it is optimal for (untransposed) matrix-vector multiplications and can be used for very large sparse matrices (not that it matters in this small example):

```csharp
// Create the matrix
IMatrix A = Matrix.CreateFromArray(new double[,] 
{
    { 10.0, 0.0, 2.0, 0.0 },
    { 0.0, 10.0, 0.0, 2.0 },
    { 2.0, 0.0, 10.0, 0.0 },
    { 0.0, 2.0, 0.0, 10.0 }
});
DokRowMajor dokA = DokRowMajor.CreateEmpty(4, 4);
for (int i = 0; i < 4; ++i)
{
    for (int j = 0; j < 4; ++j)
    {
        if (A[i, j] != 0) dokA[i, j] = A[i, j];
    }
}
CsrMatrix csrA = dokA.BuildCsrMatrix(sortRowsCols: true);

// Calculate the preconditioner
var M = new JacobiPreconditioner(A.GetDiagonalAsArray());

// PCG settings
var builder = new PcgAlgorithm.Builder();
builder.ResidualTolerance = 1E-7;
builder.MaxIterationsProvider = new FixedMaxIterationsProvider(4); // for a 4x4 matrix no more than 4 iterations should be required. Even 4 are too many. 
PcgAlgorithm pcg = builder.Build();

// Solve the linear system
Vector b = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
Vector x = Vector.CreateZero(4);
IterativeStatistics stats = pcg.Solve(csrA, M, b, x, true, () => Vector.CreateZero(b.Length));
Debug.Assert(stats.HasConverged == true);
```

# Minimum Residual
Let us solve the linear system `A * x = b`, where:

<a href="https://www.codecogs.com/eqnedit.php?latex=A&space;=&space;\begin{bmatrix}&space;2&space;&&space;&&space;&&space;\\&space;&&space;1&space;&&space;&&space;\\&space;&&space;&&space;-1&space;&&space;\\&space;&&space;&&space;&&space;-2&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A&space;=&space;\begin{bmatrix}&space;2&space;&&space;&&space;&&space;\\&space;&&space;1&space;&&space;&&space;\\&space;&&space;&&space;-1&space;&&space;\\&space;&&space;&&space;&&space;-2&space;\end{bmatrix}" title="A = \begin{bmatrix} 2 & & & \\ & 1 & & \\ & & -1 & \\ & & & -2 \end{bmatrix}" /></a>  

A is symmetric, invertible, but not positive definite. We will solve the linear system using the Minimum Residual method. The Jacobi preconditioner will be used, due to its simplicity. The matrix will be in CSR format, since it is optimal for (untransposed) matrix-vector multiplications and can be used for very large sparse matrices (not that it matters in this small example):

```csharp
// Create the matrix
IMatrix A = Matrix.CreateFromArray(new double[,] 
{
    { 2.0, 0.0,  0.0,  0.0 },
    { 0.0, 1.0,  0.0,  0.0 },
    { 0.0, 0.0, -1.0,  0.0 },
    { 0.0, 0.0,  0.0,  -2.0 }
});
DokRowMajor dokA = DokRowMajor.CreateEmpty(4, 4);
for (int i = 0; i < 4; ++i)
{
    for (int j = 0; j < 4; ++j)
    {
        if (A[i, j] != 0) dokA[i, j] = A[i, j];
    }
}
CsrMatrix csrA = dokA.BuildCsrMatrix(sortRowsCols: true);

// Calculate the preconditioner
var positiveDiagonal = new double[n];
for (int i = 0; i < n; ++i) positiveDiagonal[i] = Math.Abs(A[i, i]);
var M = new JacobiPreconditioner(positiveDiagonal);

// MINRES settings
double residualTolerance = 1E-7;
int numStoredOrthogonalDirections = 0;
bool checkMatrixSymmetry = false;
bool printIterations = false;
var minres = new MinRes(A.NumRows, residualTolerance, numStoredOrthogonalDirections, checkMatrixSymmetry, printIterations);

// Solve the linear system
Vector b = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
(IVector x, MinresStatistics stats) = minres.Solve(A, b);
Debug.Assert(stats.HasConverged == true);
```

# Generalized Minimum Residual
Let us solve the linear system `A * x = b`, where:

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\begin{bmatrix}&space;3&space;&&space;8&space;&&space;0&space;&&space;15\\&space;0&space;&&space;1&space;&&space;3&space;&&space;-7\\&space;6&space;&&space;0&space;&&space;12&space;&&space;0\\&space;3&space;&&space;0&space;&&space;9&space;&&space;15&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\begin{bmatrix}&space;3&space;&&space;8&space;&&space;0&space;&&space;15\\&space;0&space;&&space;1&space;&&space;3&space;&&space;-7\\&space;6&space;&&space;0&space;&&space;12&space;&&space;0\\&space;3&space;&&space;0&space;&&space;9&space;&&space;15&space;\end{bmatrix}" title="A=\begin{bmatrix} 3 & 8 & 0 & 15\\ 0 & 1 & 3 & -7\\ 6 & 0 & 12 & 0\\ 3 & 0 & 9 & 15 \end{bmatrix}" /></a>

A is invertible, but not symmetric. We will solve the linear system using the Generalized Minimum Residual method. No preconditioning will be used. The matrix will be in CSR format, since it is optimal for (untransposed) matrix-vector multiplications and can be used for very large sparse matrices (not that it matters in this small example):

```csharp
// Create the matrix
IMatrix A = Matrix.CreateFromArray(new double[,] 
{
    { 3.0, 8.0,  0.0,  15.0 },
    { 0.0, 1.0,  3.0,  -7.0 },
    { 6.0, 0.0, 12.0,   0.0 },
    { 3.0, 0.0,  9.0,  15.0 }
});
DokRowMajor dokA = DokRowMajor.CreateEmpty(4, 4);
for (int i = 0; i < 4; ++i)
{
    for (int j = 0; j < 4; ++j)
    {
        if (A[i, j] != 0) dokA[i, j] = A[i, j];
    }
}
CsrMatrix csrA = dokA.BuildCsrMatrix(sortRowsCols: true);

// No preconditioning
var M = new IdentityPreconditioner();

// GMRES settings
var builder = new GmresAlgorithm.Builder();
builder.MaximumIterations = 20;
builder.InnerIterationsProvider = new FixedMaxIterationsProvider(4);
builder.AbsoluteTolerance = 1E-7;
builder.RelativeTolerance = 1E-7;
GmresAlgorithm gmres = gmresAlgorithmBuilder.Build();

// Solve the linear system
Vector b = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0 });
Vector x = Vector.CreateZero(4);
gmres.Solve(matrix, new IdentityPreconditioner(), b, x, initialGuessIsZero: true, () => Vector.CreateZero(4));
Debug.Assert(stats.HasConverged == true);
```