using System;
using MGroup.Environments;

namespace MGroup.LinearAlgebra.Distributed.Tests
{
	public static class Utilities
	{
		public static IComputeEnvironment CreateEnvironment(this EnvironmentChoice environmentChoice)
		{
			if (environmentChoice == EnvironmentChoice.SequentialSharedEnvironment)
			{
				return new SequentialSharedEnvironment();
			}
			else if (environmentChoice == EnvironmentChoice.TplSharedEnvironment)
			{
				return new TplSharedEnvironment();
			}
			else
			{
				throw new NotImplementedException();
			}
		}
	}
}
