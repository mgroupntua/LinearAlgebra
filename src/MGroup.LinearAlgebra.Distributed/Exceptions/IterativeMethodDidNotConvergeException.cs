using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.LinearAlgebra.Distributed.Exceptions
{
	/// <summary>
	/// The exception that is thrown when an iterative method for solving linear systems terminates without reaching 
	/// the desired accuracy. This usually indicates that the linear system is not well defined or that the iterative method 
	/// chosen is not appropriate for that linear system. However, it is also possible that the desired accuracy is just too 
	/// strict, in which case relaxing it will allow the iterative method to converge to an acceptable solution.
	/// </summary>
	public class IterativeMethodDidNotConvergeException : Exception
	{
		/// <summary>
		/// Initializes a new instance of the <see cref="IterativeMethodDidNotConvergeException"/> class.
		/// </summary>
		public IterativeMethodDidNotConvergeException()
		{ }

		/// <summary>
		/// Initializes a new instance of the <see cref="IterativeMethodDidNotConvergeException"/> class with a specified error 
		/// message.
		/// </summary>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		public IterativeMethodDidNotConvergeException(string message) : base(message)
		{ }

		/// <summary>
		/// Initializes a new instance of the <see cref="IterativeMethodDidNotConvergeException"/> class with a specified error 
		/// message and a reference to the inner exception that is the cause of this exception.
		/// </summary>
		/// <param name="message">The error message that explains the reason for the exception.</param>
		/// <param name="inner">The exception that is the cause of the current exception. If the innerException parameter is not 
		///     a null reference, the current exception is raised in a catch block that handles the inner exception. </param>
		public IterativeMethodDidNotConvergeException(string message, Exception inner) : base(message, inner)
		{ }
	}
}
