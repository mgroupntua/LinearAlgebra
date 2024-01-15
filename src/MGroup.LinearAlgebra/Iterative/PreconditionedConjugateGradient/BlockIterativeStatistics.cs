using System;
using System.Collections.Generic;
using System.Text;

using MGroup.LinearAlgebra.Vectors;

//TODO: Add time measurements, flop measurements, etc. 
//TODO: Alternatively I could use a logger (pull observer) in the algorithm
//TODO: Needs a better name
//TODO: Each algorithm/author outputs something different. Once enough have been implemented/ported, find an appropriate design
//      to unify them. 
namespace MGroup.LinearAlgebra.Iterative
{
	/// <summary>
	/// Data Transfer Object that collects information about the execution and convergence when solving a linear system with 
	/// an iterative algorithm.
	/// </summary>
	public class BlockIterativeStatistics : IterativeStatistics
	{
		/// The square residual vector r^2
		public double ResDotRes { get; set; }

		///  The solution of the algorithm
		public IVector Solution { get; set; }
	}
}
