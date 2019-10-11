# Create vectors
Create dense vectors
```csharp
Vector x = Vector.CreateZero(10); // creates a vector with 10 zero entries
Vector y = Vector.CreateFromArray(new double[] { 10.0, 0.0, 3.0, 5.0, 6.0 }); // creates a vector with the specified array
 ```

Create sparse vectors to represent <a href="https://www.codecogs.com/eqnedit.php?latex=\begin{bmatrix}3&space;&&space;0&space;&&space;5&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;4&space;&&space;0&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix}3&space;&&space;0&space;&&space;5&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;4&space;&&space;0&space;\end{bmatrix}" title="\begin{bmatrix}3 & 0 & 5 & 0 & 0 & 0 & 4 & 0 \end{bmatrix}" /></a>
```csharp
Dictionary entries = new Dictionary();
entries[0] = 3.0;
entries[2] = 5.0;
entries[6] = 4.0;
SparseVector x = SparseVector.CreateFromDictionary(8, entries);

SparseVector y = SparseVector.CreateFromDense(new double[] { 3.0, 0.0, 5.0, 0.0, 0.0, 0.0, 4.0, 0.0 });
 ```

# Indexing
We can find the length of a vector and get or set the entry at some index:
```csharp
Vector x = Vector.CreateFromArray(new double[] { 10.0, 0.0, 3.0, 5.0, 6.0 });
SparseVector y = SparseVector.CreateFromDense(new double[] { 3.0, 0.0, 5.0, 0.0, 0.0, 0.0, 4.0, 0.0 });

int length = y.Length; // length = number of entries: y.Length = 8
double x1 = x[1]; // get entry at index
x[3] = 3.1; // set entry at index
double y1 = y[1]; // 0 will be reuturned
//y[3] = 3.0; // this is not allowed, to prevent modification of zero entries
 ```

# Common operations
Linear combinations:
```csharp
Operation                  Code
                           IVector x, y, z; // Initialize them with the same dimensions
z = x + y                  z = x.Add(y);
z = x - y                  z = x.Subtract(y);
z = 2 * x                  z = x.Scale(2);
z = y + 2 * x              z = y.Axpy(x, 2)
z = 2 * x + 3 * y          z = y.LinearCombination(3, x, 2);
```
Variations of the above can be used to overwrite one of the operands, instead of allocating new vectors. However these will fail, if they try to overwrite entries that are not explicitly stored:
```csharp
Operation                  Code
                           IVector x, y, z; // Initialize them with the same dimensions
y = x + y				   y.AddIntoThis(x);
y = x - y                  y.SubtractIntoThis(x);
x = 2 * x                  x.ScaleIntoThis(2);
y = y + 2 * x              y.AxpyIntoThis(x, 2);
y = 2 * x + 3 * y          y.LinearCombinationIntoThis(3, x, 2);
```

Dot (inner) vector product and euclidian norm:
```csharp
Operation                  Code
                           IVector x, y; // Initialize them with the same dimensions
a = x * y				   double a = x.DotProduct(y);
b = ||x||                  double b = x.Norm2();
```

# Arbitrary entrywise operations
User can choose which operation to apply to each entry:
```csharp
IVector x, y;
IVector xSquared = x.DoToAllEntries(xi => xi * xi); // Single vector: xSquared[i] = x[i] * x[i]
IVector c = x.DoEntrywise(y, (xi, yi) => (yi - xi) / Math.Min(xi, 0.1)); // Between 2 vectors: c[i] = (y[i] - x[i]) / min(x[i], 0.1)
```

# Reductions
Process all entries of a vector and calculate a single value from them. E.g. the euclidian norm can be written as follows, albeit that would be less efficient than calling `IVector.Norm()`:
```csharp
IVector x;
double identity = 0.0;
var processEntry = (entry, sum) => entry + sum;
var processZero = (numZeros, sum) => sum;
var finalize = (sum) => Math.Sqrt(sum);
x.Reduce(identity, processEntry, processZero, finalize)
```

Common reductions are available for easy use:
```csharp
IVector x;
double a = x.Sum(); // Sum of entries
double b = x.Product() // Product of entries
double c = x.Average();
double d = x.Min(); // Minimum entry
double e = x.Max();
double f = x.MinAbsolute(); // Entry that has the minimum absolute value
double g = x.MaxAbsolute();
```

# Operations at only some indices
Reading select entries:
```csharp
IVector x = Vector.CreateFromArray(new double[] { 10, 11, 12, 13, 14, 15, 16 });
Vector tail = x.GetSubvector(4, x.Length); // tail = { 14, 15, 16 }
Vector evens = x.GetSubvector(new int[]{ 0, 2, 4, 6 }); // evens = { 10, 12, 14, 16 }
```

Setting select entries:
```csharp
IVector x1 = Vector.CreateFromArray(new double[] { 10, 11, 12, 13, 14, 15, 16 });
IVector x2 = Vector.CreateFromArray(new double[] { 10, 11, 12, 13, 14, 15, 16 });
IVector x3 = Vector.CreateFromArray(new double[] { 10, 11, 12, 13, 14, 15, 16 });
IVector y = Vector.CreateFromArray(new double[] { 20, 21, 22, 23, 24, 25, 26 });
x1.CopyFrom(y); // Copies the whole vector
x2.CopyNonContiguouslyFrom(new int[] { 0, 4 }, y, new int[] { 1, 3 }); // x2 = { 21, 11, 12, 13, 23, 15, 16 }
x3.CopySubvectorFrom(4, y, 2, 3); // Copy 3 entries of y starting at index=2 to x starting at index=4: x3 = {10, 11, 12, 13, 22, 23, 24}
```

Linear operations at select entries:
```csharp
IVector x1 = Vector.CreateFromArray(new double[] { 10, 11, 12, 13, 14, 15, 16 });
IVector x2 = Vector.CreateFromArray(new double[] { 10, 11, 12, 13, 14, 15, 16 });
IVector x3 = Vector.CreateFromArray(new double[] { 10, 11, 12, 13, 14, 15, 16 });
IVector y = Vector.CreateFromArray(new double[] { 20, 21, 22, 23, 24});

x1.AddIntoThisNonContiguouslyFrom(new int[] { 0, 4 }, y, new int[] { 1, 3 }); // x1 = { 31, 11, 12, 13, 37, 15, 16 }

// Linear combination of 3 entries of y starting at index=2 to x starting at index=4: 
x2.AddSubvectorIntoThis(4, y, 2, 3); // x2 = {10, 11, 12, 13, 36, 38, 40}
x3.AxpySubvectorIntoThis(4, y, 2.0, 2, 3) // x3 = {10, 11, 12, 13, 58, 61, 64}
```

# Other operations
Join vectors:
```csharp
Vector x = Vector.CreateFromArray(new double[] { 0, 1, 2, 3 });
Vector y = Vector.CreateFromArray(new double[] { 4, 5, 6 });
Vector z = x.Append(y); // c = { 0, 1, 2, 3, 4, 5, 6 }
```

Calculate the tensor product of <a href="https://www.codecogs.com/eqnedit.php?latex=x&space;=&space;\begin{bmatrix}x_0&space;&&space;x_1\end{bmatrix}^T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x&space;=&space;\begin{bmatrix}x_0&space;&&space;x_1\end{bmatrix}^T" title="x = \begin{bmatrix}x_0 & x_1\end{bmatrix}^T" /></a>
and <a href="https://www.codecogs.com/eqnedit.php?latex=y&space;=&space;\begin{bmatrix}y_0&space;&&space;y_1&space;&&space;y_2\end{bmatrix}^T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y&space;=&space;\begin{bmatrix}y_0&space;&&space;y_1&space;&&space;y_2\end{bmatrix}^T" title="y = \begin{bmatrix}y_0 & y_1 & y_2\end{bmatrix}^T" /></a>, which is defined as

<a href="https://www.codecogs.com/eqnedit.php?latex=A&space;=&space;x&space;*&space;y^T&space;=&space;\begin{bmatrix}x_0*y_0&space;&&space;x_0*y_1&space;&&space;x_0*y_2&space;\\&space;x_1*y_0&space;&&space;x_1*y_1&space;&&space;x_1*y_2&space;\end{bmatrix}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A&space;=&space;x&space;*&space;y^T&space;=&space;\begin{bmatrix}x_0*y_0&space;&&space;x_0*y_1&space;&&space;x_0*y_2&space;\\&space;x_1*y_0&space;&&space;x_1*y_1&space;&&space;x_1*y_2&space;\end{bmatrix}" title="A = x * y^T = \begin{bmatrix}x_0*y_0 & x_0*y_1 & x_0*y_2 \\ x_1*y_0 & x_1*y_1 & x_1*y_2 \end{bmatrix}" /></a>

```csharp
Vector x = Vector.CreateFromArray(new double[] { 1, 2,});
Vector y = Vector.CreateFromArray(new double[] { 3, 4, 5 });
Matrix A = x.TensorProduct(y);
```