namespace MGroup.LinearAlgebra.Tests
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	public class NativeLibsToTest
	{
		public NativeLibsToTest(bool win64IntelMkl, bool win64SuiteSparse)
		{
			this.Win64IntelMkl = win64IntelMkl;
			this.Win64SuiteSparse = win64SuiteSparse;
		}

		public bool Win64IntelMkl { get; set; }

		public bool Win64SuiteSparse { get; set; }

		public static NativeLibsToTest CreateWithNone() => new NativeLibsToTest(false, false);
	}
}
