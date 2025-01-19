using System;
using MGroup.LinearAlgebra.Implementations;
using MGroup.LinearAlgebra.Implementations.Managed;

//TODO: These should be thread-safe. Update the documentation as well.
//TODO: A different approach is to employ the Abstract Factory pattern. Clients would use factories to create matrices, instead
//      of constructors or static factory methods. These factories would have different implementations: MKL, Managed, CUDA, etc.
//      Advantage: the client could use simultaneously more than one providers. Disadvantages: 1) A lot of methods would need 
//      references to the factory objects (this could be circumvented with singletons/ enum classes) 2) What happens for matrices
//      or vectors that need e.g. both a BLAS and a SuiteSparse provider?
//TODO: The initialization should be done at the beginning of the program to avoid interference with timing.
//TODO: The name dll in the comments of this class and the providers is too windows oriented. Once I have installed the native 
//      libraries for other operating systems, I need to stop calling them dlls.
namespace MGroup.LinearAlgebra
{
	/// <summary>
	/// Allows the user to set global settings, such as whether to use managed C# linear algebra libraries or optimized native
	/// ones. Methods exposed here are not thread-safe yet.
	/// </summary>
	public static class LibrarySettings
	{
		static LibrarySettings()
		{
			GlobalProvider = new ManagedSequentialImplementationProvider();
		}

		/// <summary>
		/// Linear algebra providers are classes that wrap calls to custom or 3rd party libraries for linear algebra operations.
		/// This property allows the user to choose which library will be used. E.g. For good performance MKL providers (native
		/// dlls) should be chosen, but if the user hasn't installed Intel MKL, then they can always fall back into the managed
		/// C# providers. The default behaviour is to use managed providers, which are comparatively inefficient. It is strongly
		/// recommended to use other providers, if the required libraries are indeed installed.
		/// </summary>
		public static IImplementationProvider GlobalProvider { get; set; }

		public static bool ThrowExceptionOnKnownPerformanceBottlenecksInReleaseBuilds { get; set; } = true;
	}
}
