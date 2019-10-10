In the examples shown on this page the following matrices will be used. 

<a href="https://www.codecogs.com/eqnedit.php?latex=A&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;\end{bmatrix}" title="A = \begin{bmatrix} 9 & 0 & 3 & 0 \\ 0 & 8 & 0 & 0 \\ 0 & 2 & 6 & 0 \\ 1 & 0 & 0 & 5 \end{bmatrix}" /></a>, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\:&space;B&space;=&space;\begin{bmatrix}&space;1&space;&&space;5&space;&&space;9&space;&&space;13&space;\\&space;2&space;&&space;6&space;&&space;10&space;&&space;14&space;\\&space;3&space;&&space;7&space;&&space;11&space;&&space;15&space;\\&space;4&space;&&space;8&space;&&space;12&space;&&space;16&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\:&space;B&space;=&space;\begin{bmatrix}&space;1&space;&&space;5&space;&&space;9&space;&&space;13&space;\\&space;2&space;&&space;6&space;&&space;10&space;&&space;14&space;\\&space;3&space;&&space;7&space;&&space;11&space;&&space;15&space;\\&space;4&space;&&space;8&space;&&space;12&space;&&space;16&space;\end{bmatrix}" title="\: B = \begin{bmatrix} 1 & 5 & 9 & 13 \\ 2 & 6 & 10 & 14 \\ 3 & 7 & 11 & 15 \\ 4 & 8 & 12 & 16 \end{bmatrix}" /></a>

Their initialization for some matrix formats is:
```csharp
Matrix A = Matrix.CreateFromArray(new double[,] 
{
    { 9.0, 0.0, 3.0, 0.0 },
    { 0.0, 8.0, 0.0, 0.0 },
    { 0.0, 2.0, 6.0, 0.0 },
    { 1.0, 0.0, 0.0, 5.0 }
});

Matrix B = Matrix.CreateFromArray(new double[,] 
{
    { 1.0, 5.0,  9.0, 13.0 },
    { 2.0, 6.0, 10.0, 14.0 },
    { 3.0, 7.0, 11.0, 15.0 },
    { 4.0, 8.0, 12.0, 16.0 }
});

double[] csrValuesA = { 9, 3, 8, 2, 6, 1, 5 };
int[] csrColIndicesA = { 0, 2, 1, 1, 2, 0, 3 };
int[] csrRowOffsetsA = { 0, 2, 3, 5, 7 };
CsrMatrix csrA = CsrMatrix.CreateFromArrays(4, 4, csrValues, csrColIndices, csrRowOffsets, true);

double[] cscValuesA = { 9, 1, 8, 2, 3, 6, 5 };
int[] cscRowIndicesA = { 0, 3, 1, 2, 0 };
int[] cscColOffsetsA = { 0, 2, 4, 6, 7 };
CscMatrix cscA = CsrMatrix.CreateFromArrays(4, 4, cscValues, cscRowIndicesA, cscColOffsetsA, true);
```

Also make sure that the matrices and vectors have the correct dimensions before applying some operation to them. In *debug* configurations that will be checked by the LinearAlgebra library as well.

# Indexing
We can find the dimensions of a matrix and get or set the entry at some row and column index:
```csharp
int m = A.NumRows;
int n = csrA.NumColumns;
double A11 = A[1, 1]; // get entry at (row, col)
A[1, 2] = 3.1; // set entry at (row, col)
double A13 = csrA[1, 3]; // 0 will be reuturned
//csrA[1, 3] = 3.0; // this is not allowed, to prevent modification of zero entries
 ```

# Entrywise operations
Linear combinations:
```csharp
Operation                  Code
                           IMatrix C;
C = A + B                  C = A.Add(B);
C = A - B                  C = A.Subtract(B);
C = 2 * A                  C = A.Scale(2);
C = B + 2 * A              C = B.Axpy(A, 2)
C = 2 * A + 3 * B          C = B.LinearCombination(3, A, 2);
```
Variations of the above can be used to overwrite one of the operands, instead of allocating new matrices. However these will fail if they try to overwrite entries of sparse matrix formats that are not explicitly stored. Therefore they should be used only if a dense matrix will be overwritten or if the 2 sparse matrices have the same sparsity pattern.
```csharp
Operation                  Code
B = A + B                  B.AddIntoThis(A);
B = A - B                  B.SubtractIntoThis(A);
A = 2 * A                  B.ScaleIntoThis(2);
B = B + 2 * A              B.AxpyIntoThis(A, 2);
B = 2 * A + 3 * B          B.LinearCombinationIntoThis(3, A, 2);
```

User can choose which operation to apply to each entry:
```csharp
IVector ASquared = A.DoToAllEntries(Aij => Aij * Aij); // Single matrix: ASquared[i, j] = A[i, j] * A[i, j]
IVector C = A.DoEntrywise(B, (Aij, Bij) => (Bij - Aij) / Math.Min(Aij, 0.1)); // Between 2 matrices: C[i, j] = (B[i, j] - A[i, j]) / min(A[i, j], 0.1)
```

# Reductions
Identical to [vector reductions](vectors.md#Reductions). Keep in mind that all entries of the matrix are processed, without taking into account which entries belong to which rows or columns. E.g. you can find the minimum entry of the whole matrix, but not the minimum entry of each column.

# Multiplications
Matrix-vector multiplications for almost all matrix formats are supported. When calling these methods, determine whether the transpose of the matrix should be used or not (default). E.g. :
```csharp
Operation                  Code
                           IVector x = Vector.CreateFromArray(new double[] { 1.1, 2.2, 3.3, 4.4 });
                           IVector y;
y = A * x                  y = A.Multiply(x);
                           y = csrA.Multiply(x); // CSR is optimal for untransposed multiplications
                           y = cscA.Multiply(x);
y = A^T * x                y = A.Multiply(x, true);
                           y = csrA.Multiply(x, true); 
                           y = cscA.Multiply(x, true); // CSC is optimal for transposed multiplications
y = x^T * A                // same as y = A^T * x
 ```

Matrix-matrix multiplications for most matrix formats are supported. When calling these methods, determine whether the transpose of the matrices should be used or not (default). E.g. :
```csharp
Operation                  Code
                           Matrix C;
C = A * B                  C = A.MultiplyRight(B);
                           C = csrA.MultiplyRight(A); // CSR on the left is optimal for untransposed multiplications
C = B * A                  C = A.MultiplyLeft(B);
                           C = cscA.MultiplyLeft(B); // CSC on the right is optimal for untransposed multiplications
C = A^T * B                C = A.MultiplyRight(B, true, false);
                           C = cscA.MultiplyRight(B, true, false); // CSC on the left is optimal for transposed multiplications
C = A * B^T                C = A.MultiplyRight(B, false, true);
C = A^T * B^T              C = A.MultiplyRight(B, true, true);
 ```

# Transposition
Usually transposing a matrix is needed during multiplication, in which case it is done implicitly by using the correct flag in the multiplication method call, as shown in [Multiplications](#Multiplications). Nevertheless, it is possible to explicitly transpose a matrix, at the cost of extra memory:
```csharp
Matrix transposeA = A.Transpose();

// In this case, no extra memory is needed. However modifying one of csrA, cscTransposeA will also modify the other.
bool copyInternalArrays = false;
CscMatrix cscTransposeA = csrA.TransposeToCSC(copyInternalArrays); 
 ```

# Join matrices
We can join 2 matrices by putting: 
- one on the left and the other on the right 

<a href="https://www.codecogs.com/eqnedit.php?latex=D&space;=&space;\begin{bmatrix}A&space;&&space;B\end{bmatrix}&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;&&space;1&space;&&space;5&space;&&space;9&space;&&space;13&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;&&space;2&space;&&space;6&space;&&space;10&space;&&space;14&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;&&space;3&space;&&space;7&space;&&space;11&space;&&space;15&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;&&space;4&space;&&space;8&space;&&space;12&space;&&space;16&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D&space;=&space;\begin{bmatrix}A&space;&&space;B\end{bmatrix}&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;&&space;1&space;&&space;5&space;&&space;9&space;&&space;13&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;&&space;2&space;&&space;6&space;&&space;10&space;&&space;14&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;&&space;3&space;&&space;7&space;&&space;11&space;&&space;15&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;&&space;4&space;&&space;8&space;&&space;12&space;&&space;16&space;\end{bmatrix}" title="D = \begin{bmatrix}A & B\end{bmatrix} = \begin{bmatrix} 9 & 0 & 3 & 0 & 1 & 5 & 9 & 13 \\ 0 & 8 & 0 & 0 & 2 & 6 & 10 & 14 \\ 0 & 2 & 6 & 0 & 3 & 7 & 11 & 15 \\ 1 & 0 & 0 & 5 & 4 & 8 & 12 & 16 \end{bmatrix}" /></a>

- one on the top and the other on the bottom

<a href="https://www.codecogs.com/eqnedit.php?latex=E&space;=&space;\begin{bmatrix}A&space;\\&space;B\end{bmatrix}&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;\\&space;1&space;&&space;5&space;&&space;9&space;&&space;13&space;\\&space;2&space;&&space;6&space;&&space;10&space;&&space;14&space;\\&space;3&space;&&space;7&space;&&space;11&space;&&space;15&space;\\&space;4&space;&&space;8&space;&&space;12&space;&&space;16&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?E&space;=&space;\begin{bmatrix}A&space;\\&space;B\end{bmatrix}&space;=&space;\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;&&space;0&space;\\&space;0&space;&&space;8&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;&&space;0&space;\\&space;1&space;&&space;0&space;&&space;0&space;&&space;5&space;\\&space;1&space;&&space;5&space;&&space;9&space;&&space;13&space;\\&space;2&space;&&space;6&space;&&space;10&space;&&space;14&space;\\&space;3&space;&&space;7&space;&&space;11&space;&&space;15&space;\\&space;4&space;&&space;8&space;&&space;12&space;&&space;16&space;\end{bmatrix}" title="E = \begin{bmatrix}A \\ B\end{bmatrix} = \begin{bmatrix} 9 & 0 & 3 & 0 \\ 0 & 8 & 0 & 0 \\ 0 & 2 & 6 & 0 \\ 1 & 0 & 0 & 5 \\ 1 & 5 & 9 & 13 \\ 2 & 6 & 10 & 14 \\ 3 & 7 & 11 & 15 \\ 4 & 8 & 12 & 16 \end{bmatrix}" /></a>


```csharp
Matrix D = A.AppendRight(B);
Marix E = A.AppendBottom(B);
 ```

# Operations at only some rows & columns
Extract a single row, column or the diagonal (only for square matrices):
```csharp
Vector row1 = A.GetRow(1);
Vector col2 = A.GetColumn(2);
Vector diagonal = A.GetDiagonal();
 ```

 Modify a row or column or part of it:
 ```csharp
// Set the whole row 1
int row1 = 1;
int start = 0;
Vector newValues = Vector.CreateFromArray(new double[] { 1.1, 2.2, 3.3, 4.4 };)
A.SetSubrow(row1, newValues, start);

// Set the second half of column 2
int col2 = 2;
start = 2;
newValues = Vector.CreateFromArray(new double[] { 3.3, 4.4 };)
A.SetSubcolumn(col2, newValues, start);
 ```

 Extract the submatrices

 <a href="https://www.codecogs.com/eqnedit.php?latex=A_{sub1}&space;=&space;\begin{bmatrix}&space;A_{00}&space;&&space;A_{01}&space;&&space;A_{02}&space;\\&space;A_{10}&space;&&space;A_{11}&space;&&space;A_{12}&space;\\&space;A_{20}&space;&&space;A_{21}&space;&&space;A_{22}&space;\\&space;\end{bmatrix}&space;=\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;\\&space;0&space;&&space;8&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;\\&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A_{sub1}&space;=&space;\begin{bmatrix}&space;A_{00}&space;&&space;A_{01}&space;&&space;A_{02}&space;\\&space;A_{10}&space;&&space;A_{11}&space;&&space;A_{12}&space;\\&space;A_{20}&space;&&space;A_{21}&space;&&space;A_{22}&space;\\&space;\end{bmatrix}&space;=\begin{bmatrix}&space;9&space;&&space;0&space;&&space;3&space;\\&space;0&space;&&space;8&space;&&space;0&space;\\&space;0&space;&&space;2&space;&&space;6&space;\\&space;\end{bmatrix}" title="A_{sub1} = \begin{bmatrix} A_{00} & A_{01} & A_{02} \\ A_{10} & A_{11} & A_{12} \\ A_{20} & A_{21} & A_{22} \\ \end{bmatrix} =\begin{bmatrix} 9 & 0 & 3 \\ 0 & 8 & 0 \\ 0 & 2 & 6 \\ \end{bmatrix}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=A_{sub2}&space;=&space;\begin{bmatrix}&space;A_{33}&space;&&space;A_{30}&space;\\&space;A_{03}&space;&&space;A_{00}&space;\\&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;5&space;&&space;1&space;\\&space;0&space;&&space;9&space;\\&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A_{sub2}&space;=&space;\begin{bmatrix}&space;A_{33}&space;&&space;A_{30}&space;\\&space;A_{03}&space;&&space;A_{00}&space;\\&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;5&space;&&space;1&space;\\&space;0&space;&&space;9&space;\\&space;\end{bmatrix}" title="A_{sub2} = \begin{bmatrix} A_{33} & A_{30} \\ A_{03} & A_{00} \\ \end{bmatrix} = \begin{bmatrix} 5 & 1 \\ 0 & 9 \\ \end{bmatrix}" /></a>

 ```csharp
// Extract submatrix from (0, 0) to (1, 2)
Matrix Asub1 = A.GetSubmatrix(0, 2, 0, 3);

// Extract rows and columns 0, 3, but also reverse them:
Matrix Asub2 = A.GetSubmatrix(new int[] { 3, 0 }, new int[] { 3, 0 });
 ```

Set all entries of the submatrix A<sub>sub1</sub> to 3.14:

```csharp
Matrix newValues = Matrix.CreateZero(3, 3);
newValues.SetAll(3.14);
Matrix Asub1 = A.SetSubmatrix(0, 0, newValues);
 ```
