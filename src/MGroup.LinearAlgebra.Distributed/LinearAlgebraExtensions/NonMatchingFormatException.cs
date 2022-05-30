using System;
using System.Collections.Generic;
using System.Text;

//TODO: This belongs in LinearAlgebra project
namespace MGroup.LinearAlgebra.Distributed.LinearAlgebraExtensions
{
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
