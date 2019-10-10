# Dense matrices
When working with dense matrix formats, such as e.g. Matrix, SymmetricMatrix, TriangularLower, TriangularUpper, creating them is straightforward. E.g. create the following matrices:

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\begin{bmatrix}&space;1&space;&&space;4&space;&&space;7\\&space;2&space;&&space;5&space;&&space;8\\&space;3&space;&&space;6&space;&&space;9&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\begin{bmatrix}&space;1&space;&&space;4&space;&&space;7\\&space;2&space;&&space;5&space;&&space;8\\&space;3&space;&&space;6&space;&&space;9&space;\end{bmatrix}" title="A=\begin{bmatrix} 1 & 4 & 7\\ 2 & 5 & 8\\ 3 & 6 & 9 \end{bmatrix}" /></a> , 
<a href="https://www.codecogs.com/eqnedit.php?latex=B=\begin{bmatrix}&space;1&space;&&space;2&space;&&space;4\\&space;2&space;&&space;3&space;&&space;5\\&space;4&space;&&space;5&space;&&space;6&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?B=\begin{bmatrix}&space;1&space;&&space;2&space;&&space;4\\&space;2&space;&&space;3&space;&&space;5\\&space;4&space;&&space;5&space;&&space;6&space;\end{bmatrix}" title="B=\begin{bmatrix} 1 & 2 & 4\\ 2 & 3 & 5\\ 4 & 5 & 6 \end{bmatrix}" /></a> , 
<a href="https://www.codecogs.com/eqnedit.php?latex=C=\begin{bmatrix}&space;1&space;&&space;2&space;&&space;4\\&space;0&space;&&space;3&space;&&space;5\\&space;0&space;&&space;0&space;&&space;6&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?C=\begin{bmatrix}&space;1&space;&&space;2&space;&&space;4\\&space;0&space;&&space;3&space;&&space;5\\&space;0&space;&&space;0&space;&&space;6&space;\end{bmatrix}" title="C=\begin{bmatrix} 1 & 2 & 4\\ 0 & 3 & 5\\ 0 & 0 & 6 \end{bmatrix}" /></a>

```csharp
Matrix A = Matrix.CreateZero(3, 3);
double temp = 0;
for (int j = 0; j < 3; ++j)
{
    for (int i = 0; i < 3; ++i)
    {
        ++temp;
        A[i, j] = temp;
    }
}

SymmetricMatrix B = SymmetricMatrix.CreateZero(3);
TriangularUpper C = TriangularUpper.CreateZero(3);
temp = 0;
for (int j = 0; j < 3; ++j)
{
    for (int i = 0; i <= j; ++i)
    {
        ++temp;
        B[i, j] = temp;
        C[i, j] = temp;
    }
}
```

# Sparse matrices
Sparse matrix formats that are efficient for certain linear algebra operations (e.g. CSR for multiplication) are unsuitable for randomly setting entries. Their sparsity pattern must be determined before they are created. For this reason DOK sparse matrix formats, which support efficient setting of random entries, are provided. After creating the matrix in those, convert them into the sparse matrix formats needed in each application. E.g. create the following matrix in various formats:

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;6&space;\\&space;&&space;2&space;&&space;5&space;&&space;0&space;\\&space;&&space;&&space;3&space;&&space;0&space;\\&space;sym&space;&&space;&&space;&&space;4&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;6&space;\\&space;&&space;2&space;&&space;5&space;&&space;0&space;\\&space;&&space;&&space;3&space;&&space;0&space;\\&space;sym&space;&&space;&&space;&&space;4&space;\end{bmatrix}" title="A=\begin{bmatrix} 1 & 0 & 0 & 6 \\ & 2 & 5 & 0 \\ & & 3 & 0 \\ sym & & & 4 \end{bmatrix}" /></a>

```csharp
DokColMajor dokForCsc = DokColMajor.CreateEmpty(4, 4);
DokRowMajor dokForCsr = DokRowMajor.CreateEmpty(4, 4);
DokSymmetric dokForSymCsc = DokSymmetric.CreateEmpty(4);

// Diagonal entries
double temp = 0;
for (int i = 0; i < 4; ++i)
{
    ++temp;
    dokForCsc[i, i] = temp;
    dokForCsr[i, i] = temp;
    dokForSymmCsc[i, i] = temp;
}

// Other entries
dokForCsc[1, 2] = 5.0;
dokForCsc[2, 1] = 5.0;
dokForCsr[1, 2] = 5.0;
dokForCsr[2, 1] = 5.0;
dokForSymCsc[1, 2] = 5.0;
dokForCsc[0, 3] = 6.0;
dokForCsc[3, 0] = 6.0;
dokForCsr[0, 3] = 6.0;
dokForCsr[0, 3] = 6.0;
dokForSymCsc[0, 3] = 6.0;

// Convert DOKs into the matrix formats used for operations
bool sortRowsCols = true;
CscMatrix csc = dokForCsc.BuildCscMatrix(sortRowsCols);
CsrMatrix csr = dokForCsr.BuildCsrMatrix(sortRowsCols);
SymmetricCscMatrix symCsc = dokForSymCsc.BuildSymmetricCscMatrix(sortRowsCols);
```

# Finite element matrix assembly
In the Finite Element Method and similar partial differential equation numerical methods, we often need to assemble a global matrix from smaller (local) submatrices. E.g. assemble the global matrix A<sub>glob</sub> from three submatrices A<sub>loc</sub>:

<a href="https://www.codecogs.com/eqnedit.php?latex=A_{1\:&space;loc}=A_{2\:&space;loc}=A_{3\:&space;loc}=\begin{bmatrix}&space;2&space;&&space;1\\&space;1&space;&&space;2&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A_{1\:&space;loc}=A_{2\:&space;loc}=A_{3\:&space;loc}=\begin{bmatrix}&space;2&space;&&space;1\\&space;1&space;&&space;2&space;\end{bmatrix}" title="A_{1\: loc}=A_{2\: loc}=A_{3\: loc}=\begin{bmatrix} 2 & 1\\ 1 & 2 \end{bmatrix}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=A_{glob}=\begin{bmatrix}&space;2&space;&&space;1&space;&&space;&&space;\\&space;1&space;&&space;2&plus;{\color{Red}&space;2}&space;&&space;{\color{Red}&space;1}&space;&&space;\\&space;&&space;{\color{Red}&space;1}&space;&&space;{\color{Red}&space;2}&plus;{\color{Blue}&space;2}&space;&&space;{\color{Blue}&space;1}&space;\\&space;&&space;&&space;{\color{Blue}&space;1}&space;&&space;{\color{Blue}&space;{\color{Blue}&space;}2}&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A_{glob}=\begin{bmatrix}&space;2&space;&&space;1&space;&&space;&&space;\\&space;1&space;&&space;2&plus;{\color{Red}&space;2}&space;&&space;{\color{Red}&space;1}&space;&&space;\\&space;&&space;{\color{Red}&space;1}&space;&&space;{\color{Red}&space;2}&plus;{\color{Blue}&space;2}&space;&&space;{\color{Blue}&space;1}&space;\\&space;&&space;&&space;{\color{Blue}&space;1}&space;&&space;{\color{Blue}&space;{\color{Blue}&space;}2}&space;\end{bmatrix}" title="A_{glob}=\begin{bmatrix} 2 & 1 & & \\ 1 & 2+{\color{Red} 2} & {\color{Red} 1} & \\ & {\color{Red} 1} & {\color{Red} 2}+{\color{Blue} 2} & {\color{Blue} 1} \\ & & {\color{Blue} 1} & {\color{Blue} {\color{Blue} }2} \end{bmatrix}" /></a>

This can be done easily and efficiently using DOK matrices or other matrix builders, with the use of arrays that map the rows/columns of A<sub>loc</sub> to the rows/columns of A<sub>glob</sub>:

```csharp

// Local matrices. They are symmetric.
var localMatrices = new List<Matrix>();
for (int i = 0; i < 3; ++i)
{
    Matrix localA = Matrix.CreateZero(2, 2);
    localA[0, 0] = 2.0;
    localA[0, 1] = 1.0;
    localA[1, 0] = 1.0;
    localA[1, 1] = 2.0;
    localMatrices[i] = localA;
}

// Mapping arrays. In this case rows and columns are identical, since all matrices is symmetric.
var localIndices = new List<int[]>();
for (int i = 0; i < 3; ++i) localIndices[i] = new int[] { 0, 1 };
var globalIndices = new List<int[]>();
globalIndices[0] = new int[] { 0, 1 };
globalIndices[0] = new int[] { 1, 2 };
globalIndices[0] = new int[] { 2, 3 };

// Skyline matrices need some more data to determine the bandwidth.
int[] colHeights = { 0, 1, 1, 1 }; // Number of entries from the diagonal (exclusive) to the top non-zero entry. 

// Assemble local matrices into global
DokColMajor dokForCsc = DokColMajor.CreateEmpty(4, 4);
DokRowMajor dokForCsr = DokRowMajor.CreateEmpty(4, 4);
DokSymmetric dokForSymCsc = DokSymmetric.CreateEmpty(4);
SkylineBuilder builderForSky = SkylineBuilder.Create(4, colHeights);
for (int i = 0; i < 3; ++i)
{
    dokForCsc.AddSubmatrix(localMatrices[i], localIndices[i], globalIndices[i], localIndices[i], globalIndices[i]);
    dokForCsr.AddSubmatrix(localMatrices[i], localIndices[i], globalIndices[i], localIndices[i], globalIndices[i]);
    dokForSymCsc.AddSubmatrixSymmetric(localMatrices[i], localIndices[i], globalIndices[i]);
    builderForSky.AddSubmatrixSymmetric(localMatrices[i], localIndices[i], globalIndices[i]);
}

// Convert DOKs into the matrix formats used for operations
bool sortRowsCols = true;
CscMatrix cscGlobalA = dokForCsc.BuildCscMatrix(sortRowsCols);
CsrMatrix csrGlobalA = dokForCsr.BuildCsrMatrix(sortRowsCols);
SymmetricCscMatrix symCscGlobalA = dokForSymCsc.BuildSymmetricCscMatrix(sortRowsCols);
SkylineMatrix skyGlobalA = builderForSky.BuildSkylineMatrix();
```

