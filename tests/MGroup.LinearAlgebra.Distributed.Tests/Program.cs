using System;

using MGroup.LinearAlgebra.Distributed.Tests.Overlapping;

namespace MGroup.LinearAlgebra.Distributed.Tests
{
    public class Program
    {
        public static void Main(string[] args)
        {
            DistributedOverlappingIndexerTests.RunMpiTests();
            DistributedOverlappingVectorTests.RunMpiTests();
            DistributedOverlappingMatrixTests.RunMpiTests();
            DistributedOverlappingPcgTests.RunMpiTests();
		}
    }
}
