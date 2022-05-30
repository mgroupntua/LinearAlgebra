
![alt text](http://mgroup.ntua.gr/wp-content/uploads/2018/05/MGroup52.png "MGroup")

# LinearAlgebra
Library that provides implementations or wrappers for linear algebra operations in C#.

[![Build Status](https://dev.azure.com/mgroupntua/LinearAlgebra/_apis/build/status/mgroupntua.LinearAlgebra?branchName=develop)](https://dev.azure.com/mgroupntua/LinearAlgebra/_build/latest?definitionId=1&branchName=develop)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=LinearAlgebra&metric=alert_status)](https://sonarcloud.io/dashboard?id=LinearAlgebra)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## Features

- **Basic matrix and vector operations:** What you would expect from a BLAS library.
  * Vector-vector operations: Addition, subtraction, linear combinations, scaling, dot (inner) product
  * Matrix-vector operations: matrix-vector multiplication, back and forward substitution
  * Matrix-matrix operations: Addition, subtraction, linear combinations, scaling, matrix-matrix multiplication

- **Sparse matrix and vector formats:** Sparse formats avoid storing and processing zero entries, which reduces memory consumption and computational time.
  * Compressed Sparse Rows, Compressed Sparse Columns: Optimal for matrix-vector and matrix-matrix multiplications 
  * Dictionary Of Keys: Optimal for creating a sparse matrix and then copying it to another format.
  * Skyline: Optimal for Cholesky/LDL factorization

- **Direct linear system solvers:** Solving linear systems using factorization.
  * LU factorization with pivoting for general matrices
  * Cholesky/LDL factorization for symmetric positive definite matrices
  * Reordering algorithms to reduce the fill-in of the factorization process, thereby reducing memory consumption and computational time
  
- **Iterative linear system solvers:** Solving large sparse linear systems using algorithms based on matrix-vector multiplication and vector-vector operations.
  * Generalized Minimum Residual method for general matrices
  * Minimum Residual method for symmetric matrices
  * Conjugate Gradient method for symmetric positive definite matrices
  * General purpose preconditioners

- **Least squares problems:** Solving least squares and minimum norm problems
  * QR and LQ factorizations of rectangular matrices
  * Only for dense matrices at the moment

- **Matrix and vector I/O:** Reading matrices and vectors from or writing them to files
  * File formats for dense or sparse matrices
  * Compatible with matrix market file formats
  * Compatible with matlab file formats

- **Optimized native libraries:** Assigning operations to high performance linear algebra libraries
  * [Intel Math Kernel Library](https://software.intel.com/en-us/mkl). Mainly for dense linear algebra operations.
  * [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html). For sparse direct solvers.
  * Only for win64 operating systems at the moment.
  * These are optional. If the user does not want to install these libraries or does not use a win64 OS, there are is an alternative, less efficient C# implementation of each feature.
  


## Installation instructions
You can choose either to clone the solution or downloads it as a zip file.

### Clone solution
1. Under the repository name, click **Clone or Download** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. In the popup appearing choose the **Use HTTPS** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/2.png "2")

3. Use the ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/3.png "3") to copy the link provided.

4. Open Visual Studio. In Team Explorer window appearing in your screen under Local Git Repositories click the **Clone** option. If Team Explorer window is not visible you can enable in View -> Team Explorer

  ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/4.png "4")
  
5. In the text box appearing paste the link.

 ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/5.png "5")

6. Click clone and Visual Studio will automatically download and import **MGroup.LinearAlgebra**


### Download as ZIP
1. Under the repository name, click **Clone or Download** option

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. Click **Download ZIP** option. **MGroup.LinearAlgebra** will be downloaded as a ZIP file.

3. Extract the ZIP file to the folder of choice.

4. Double click on **MGroup.LinearAlgebra.sln** file to open the code with Visual Studio


## Reference manual
For more information on the code functionality please refer to the reference manual that presents an in depth documentation:
https://github.com/mgroupntua/LinearAlgebra/wiki
