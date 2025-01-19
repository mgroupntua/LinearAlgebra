namespace MGroup.LinearAlgebra.Exceptions
{
	using System;

	/// <summary>
	/// The exception that is thrown when a call to SuiteSparse library fails.
	/// </summary>
	public class NativeLibException : Exception
	{
		/// <summary>
		/// Initializes a new instance of the <see cref="NativeLibException"/> class with a specified error message.
		/// </summary>
		/// <param name="libName">The name of the native library.</param>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		public NativeLibException(string libName, string message)
			: base($"Native library {libName} exception: {message}")
		{ }

		/// <summary>
		/// Initializes a new instance of the <see cref="NativeLibException"/> class with a specified error message 
		/// and a reference to the inner exception that is the cause of this exception.
		/// </summary>
		/// <param name="libName">The name of the native library.</param>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		/// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
		///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
		public NativeLibException(string libName, string message, Exception inner)
			: base($"Native library {libName} exception: {message}", inner)
		{ }
	}
}
