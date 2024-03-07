using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.LinearAlgebraExtensions
{
	/// <summary>
	/// The exception that is thrown when a linear algebra operation cannot be executed due to the arguments (vectors, matrices 
	/// or both) having incompatible storage formats.
	/// </summary>
	public class NonMatchingFormatException : ArgumentException
	{
		public NonMatchingFormatException()
		{
		}

		public NonMatchingFormatException(string message)
			: base(message)
		{
		}

		public NonMatchingFormatException(string message, Exception inner)
			: base(message, inner)
		{
		}
	}
}
