using System;
using System.IO;
using MGroup.LinearAlgebra.Matrices.Operators;
using MGroup.LinearAlgebra.Output;
using MGroup.LinearAlgebra.Tests.Matrices;
using MGroup.LinearAlgebra.Tests.TestData;
using MGroup.LinearAlgebra.Tests.Utilities;
using Xunit;

namespace MGroup.LinearAlgebra.Tests.Output
{
    /// <summary>
    /// Tests for <see cref="BooleanMatrixWriter"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class BooleanMatrixWriterTests
    {
        [Fact]
        private static void TestWriter1()
        {
            var matrix = SignedBooleanMatrixTests.CreateMatrix(SignedBoolean5by10.A1);
            TestWriteOperation(matrix, SignedBoolean5by10.FilePath1);
        }

        [Fact]
        private static void TestWriter2()
        {
            var matrix = SignedBooleanMatrixTests.CreateMatrix(SignedBoolean5by10.A2);
            TestWriteOperation(matrix, SignedBoolean5by10.FilePath2);
        }

        private static void TestWriteOperation(SignedBooleanMatrixRowMajor matrix, string referenceFile)
            => TestWriteOperation(matrix, referenceFile, new BooleanMatrixWriter(true));

        private static void TestWriteOperation(SignedBooleanMatrixRowMajor matrix, string referenceFile, BooleanMatrixWriter writer)
        {
            string tempFile = Guid.NewGuid().ToString() + ".txt";
            writer.WriteToFile(matrix, tempFile);
            bool success = IOUtilities.AreFilesEquivalent(referenceFile, tempFile);
            File.Delete(tempFile);
            Assert.True(success);
        }
    }
}
