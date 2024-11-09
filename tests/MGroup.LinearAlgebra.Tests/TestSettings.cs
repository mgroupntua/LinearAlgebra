namespace MGroup.LinearAlgebra.Tests
{
	using System;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Implementations.MKL;

	using Xunit;

	// Currently SuiteSparse dlls call MKL dll
	public enum TestSuiteSparseAndMklLibs
	{
		Neither, MklOnly, Both
	}

	public class TestSettings
	{
		// Set the appropriate enums and flags here, in order to choose which native library tests will be run.
		private static readonly TestSuiteSparseAndMklLibs librariesToTest = TestSuiteSparseAndMklLibs.Neither;

		public const string MessageWhenSkippingMKL = "MKL is not set to be tested. See TestSettings.cs for more.";

		public const string MessageWhenSkippingSuiteSparse
			= "SuiteSparse is not set to be tested. See TestSettings.cs for more.";

		public static TheoryData<IImplementationProvider> ProvidersToTest
		{
			get
			{
				var theoryData = new TheoryData<IImplementationProvider>();
				theoryData.Add(new ManagedSequentialImplementationProvider());
				if ((librariesToTest == TestSuiteSparseAndMklLibs.MklOnly)
					|| (librariesToTest == TestSuiteSparseAndMklLibs.Both))
				{
					theoryData.Add(new NativeWin64ImplementationProvider());
				}

				return theoryData;
			}
		}

		public static bool TestMkl => (librariesToTest == TestSuiteSparseAndMklLibs.MklOnly)
			|| (librariesToTest == TestSuiteSparseAndMklLibs.Both);

		public static bool TestSuiteSparse => (librariesToTest == TestSuiteSparseAndMklLibs.Both);

		public static void RunMultiproviderTest(IImplementationProvider provider, Action test)
		{
			IImplementationProvider defaultProvider = LibrarySettings.GlobalProvider; // Store it for later
			LibrarySettings.GlobalProvider = provider;

			try
			{
				test();
			}
			finally
			{
				LibrarySettings.GlobalProvider = defaultProvider; // Once finished, reset the default providers
			}
		}
	}
}
