namespace MGroup.LinearAlgebra.Tests
{
	using System;
	using System.IO;
	using System.Reflection;
	using System.Text.Json;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.LinearAlgebra.Implementations.NativeWin64.MKL;
	using MGroup.LinearAlgebra.Triangulation;

	using Xunit;

	public static class TestSettings
	{
		// Explicit static constructor to tell C# compiler not to mark type as beforefieldinit. Only required for laziness.
		static TestSettings()
		{
			LibsToTest = NativeLibsToTest.CreateWithNone();
			ProvidersToTest = new TheoryData<IImplementationProvider>();
			ProvidersToTest.Add(new ManagedSequentialImplementationProvider());

			try
			{
				// Read from JSON
				string execDirectory = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);
				string jsonFile = Path.Combine(execDirectory, "NativeLibsToTest.json");
				string jsonText = File.ReadAllText(jsonFile);

				var options = new JsonSerializerOptions { PropertyNameCaseInsensitive = true };
				NativeLibsToTest deserialized = JsonSerializer.Deserialize<NativeLibsToTest>(jsonText, options);
				if (deserialized != null)
				{
					LibsToTest = deserialized;
					if (LibsToTest.Win64IntelMkl)
					{
						if (LibsToTest.Win64SuiteSparse)
						{
							ProvidersToTest.Add(new NativeWin64ImplementationProvider());
						}
						else
						{
							ProvidersToTest.Add(new CustomImplementationProvider(
								MklBlasProvider.UniqueInstance,
								MklSparseBlasProvider.UniqueInstance,
								MklLapackProvider.UniqueInstance,
								new ManagedReorderingProvider(),
								superNodal => new CholeskyCSparseNet()));
						}
					}
				}
			}
			catch (Exception ex)
			{
				// If reading native lib options fails for any reason, do nothing (use only managed providers).
			}
		}

		public static NativeLibsToTest LibsToTest { get; }

		public static TheoryData<IImplementationProvider> ProvidersToTest { get; }

		public const string SkipMessage =
			"This native library is not set to be tested. You can set it in NativeLibsToTest.json. See TestSettings.cs for more.";

		public static void RunMultiproviderTest(IImplementationProvider provider, Action test)
		{
			//TODO: The default provider logic is not thread-safe
			IImplementationProvider defaultProvider = LibrarySettings.GlobalProvider; // Store it for later.
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
