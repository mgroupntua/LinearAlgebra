namespace MGroup.LinearAlgebra.Exceptions
{
	using System;
	using System.Collections.Generic;
	using System.Text;

	/// <summary>
	/// The exception that is thrown when a matrix does not have the specific sparsity pattern, required for the desired linear 
	/// algebra operation.
	/// </summary>
	public class InvalidSparsityPatternException : Exception
	{
		/// <summary>
		/// Initializes a new instance of the <see cref="InvalidSparsityPatternException"/> class with a specified error message.
		/// </summary>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		public InvalidSparsityPatternException(string message) : base(message)
		{ }

		/// <summary>
		/// Initializes a new instance of the <see cref="InvalidSparsityPatternException"/> class with a specified error message 
		/// and a reference to the inner exception that is the cause of this exception.
		/// </summary>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		/// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
		///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
		public InvalidSparsityPatternException(string message, Exception inner) : base(message, inner)
		{ }
	}
}
