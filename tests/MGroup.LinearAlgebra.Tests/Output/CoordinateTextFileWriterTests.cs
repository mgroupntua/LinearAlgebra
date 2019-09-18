using System;
using System.IO;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Output;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="CoordinateTextFileWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class CoordinateTextFileWriterTests
    {
        [Fact]
        private static void TestSparsePosDef()
        {
            var matrix = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order, SparsePosDef10by10.SkylineValues,
                SparsePosDef10by10.SkylineDiagOffsets, true, true);
            TestWriteOperation(matrix, SparsePosDef10by10.CoordinateFormatPath);
        }

        [Fact]
        private static void TestSparseRectangularCSC()
        {
            // The reference file has been written by iterating the entries in CSC order.
            // Therefore it cannot be used for checking with a CSR matrix.
            var matrix = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            TestWriteOperation(matrix, SparseRectangular10by5.CoordinateFormatPath);
        }

        private static void TestWriteOperation(ISparseMatrix matrix, string referenceFile)
            => TestWriteOperation(matrix, referenceFile, new CoordinateTextFileWriter());

        private static void TestWriteOperation(ISparseMatrix matrix, string referenceFile, CoordinateTextFileWriter writer)
        {
            string tempFile = Guid.NewGuid().ToString() + ".txt";
            writer.WriteToFile(matrix, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
