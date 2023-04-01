using System;

namespace MGroup.LinearAlgebra.Commons
{
    public class ValueComparer
    {
        private readonly double tolerance;

        public ValueComparer(double tolerance)
        {
            if (tolerance < 0.0 || tolerance > 1.0)
                throw new ArgumentException("tolerance must belong to the interval [0,1]");
            this.tolerance = tolerance;
        }

        public bool AreEqual(double val1, double val2)
        {
            if (double.IsNaN(val1) || double.IsNaN(val2)) return false;
            if (Math.Abs(val2) <= tolerance) // Can't divide with expected ~= 0. 
            {
                if (Math.Abs(val1) <= tolerance) return true;
                else return false;
            }
            else return (Math.Abs(1.0 - val1 / val2) < tolerance) ? true : false;
        }
    }
}
