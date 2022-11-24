using System;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.Differentiation;

namespace Lab1.OptimizaciyaFunkciyOdnoyPeremennoy
{
    public class Functions
    {
        public int N { get; private set; }

        public double X { get; private set; }

        public Functions(int n = 1)
        {
            this.N = n;

            switch (this.N)
            {
                case 1:
                    this.X = -2.3247;
                    break;
                default:
                    throw new ArgumentOutOfRangeException("No such function");
            }
        }

        public double Solve(double x)
        {
            switch (this.N)
            {
                case 1:
                    return Math.Pow(x + 1, 4) - 2 * Math.Pow(x, 2);
                default:
                    throw new ArgumentOutOfRangeException("No such function");
            }
        }
    }

    public abstract class Methods
    {
        public Functions F { get; protected set; }

        public double A { get; private set; }

        public double B { get; private set; }

        public double EPS { get; private set; }

        public double X { get; protected set; }

        public double Y { get; protected set; }

        public double E { get; protected set; }

        public double I { get; protected set; }

        public double It { get; private set; }

        public double Acc { get; protected set; }

        public Methods(int n, double a, double b, double eps, int it)
        {
            F = new Functions(n);
            A = a;
            B = b;
            I = 0;
            It = it;
            if (eps <= 0 || eps >= 1)
            {
                throw new ArgumentOutOfRangeException("0 <= eps <= 1");
            }
            else
            {
                EPS = eps;
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append($"\n\tFunction №{F.N}\n " +
                $"\tMethod {GetType()}\n" +
                $"\tInterval: ({A}, {B})\n" +
                $"\tSpecified error: {EPS}\n" +
                $"\tx = {X:f10}\n" +
                $"\tf(x) = {Y:f10}\n" +
                $"\tEps = {E:f10}\n" +
                $"\tNumber of iterations: {I}\n" +
                $"\tAccuracy at {It} iterations: {Acc}");

            return sb.ToString();
        }
    }

    public class GoldenSection : Methods
    {
        public GoldenSection(int n = 1, double a = -3, double b = -2, double eps = 0.0001, int it = 4) : base(n, a, b, eps, it) 
        {
            double y = a + (3 - Math.Sqrt(5)) * (b - a) / 2;
            double z = a + b - y;

            while (Math.Abs(b - a) > EPS)
            {
                if (F.Solve(y) <= F.Solve(z))
                {
                    b = z;
                    z = y;
                    y = a + b - y;
                }
                else
                {
                    a = y;
                    y = z;
                    z = a + b - z;
                }

                I++;
            }

            X = (a + b) / 2;
            Y = F.Solve(X);
            E = Math.Abs(b - a);


            a = A;
            b = B;
            y = a + (3 - Math.Sqrt(5)) * (b - a) / 2;
            z = a + b - y;
            int i = 0;

            while (i < It)
            {
                if (F.Solve(y) <= F.Solve(z))
                {
                    b = z;
                    z = y;
                    y = a + b - y;
                }
                else
                {
                    a = y;
                    y = z;
                    z = a + b - z;
                }

                i++;
            }

            Acc = Math.Abs(b - a);
        }
    }

    public class Dichotomy : Methods
    {
        public Dichotomy(int n = 1, double a = -3, double b = -2, double eps = 0.0001, int it = 4) : base(n, a, b, eps, it)
        {
            Random rnd = new Random();
            double e = rnd.NextDouble() * eps * 2;
            double y;
            double z;

            while (Math.Abs(b - a) > EPS)
            {
                y = (a + b - e) / 2;
                z = (a + b + e) / 2;

                if (F.Solve(y) <= F.Solve(z))
                {
                    b = z;
                }
                else
                {
                    a = y;
                }

                I++;
            }

            X = (a + b) / 2;
            Y = F.Solve(X);
            E = Math.Abs(b - a);


            a = A;
            b = B;
            int i = 0;

            while (i < It)
            {
                y = (a + b - e) / 2;
                z = (a + b + e) / 2;

                if (F.Solve(y) <= F.Solve(z))
                {
                    b = z;
                }
                else
                {
                    a = y;
                }

                i++;
            }

            Acc = Math.Abs(b - a);
        }
    }

    public class Newton_Raphson : Methods
    {
        public Newton_Raphson(int n = 1, double a = -3, double b = -2, double eps = 0.0001, int it = 4) : base(n, a, b, eps, it)
        {
            Func<double, double> fx = F.Solve;
            double x1;

            if (Differentiate.Derivative(fx, A, 1) * Differentiate.Derivative(fx, A, 3) > 0)
            {
                x1 = A;
            }
            else
            {
                x1 = B;
            }

            double x2 = x1;

            while (Math.Abs(Differentiate.Derivative(fx, x1, 1)) > EPS)
            {
                x1 = x1 - Differentiate.Derivative(fx, x1, 1) / Differentiate.Derivative(fx, x1, 2);

                I++;
            }

            X = x1;
            Y = F.Solve(X);
            E = Math.Abs(Differentiate.Derivative(fx, x1, 1));


            int i = 0;

            while (i < It)
            {
                x2 = x2 - Differentiate.Derivative(fx, x2, 1) / Differentiate.Derivative(fx, x2, 2);

                i++;
            }

            Acc = Math.Abs(Differentiate.Derivative(fx, x1, 1));
        }
    }

    public class MiddlePoints : Methods
    {
        public MiddlePoints(int n = 1, double a = -3, double b = -2, double eps = 0.0001, int it = 4) : base(n, a, b, eps, it)
        {
            Func<double, double> fx = F.Solve;
            double z;

            do
            {
                z = (a + b) / 2;

                if (Differentiate.Derivative(fx, z, 1) < 0)
                {
                    a = z;
                }
                else
                {
                    b = z;
                }

                I++;
            } while (Math.Abs(Differentiate.Derivative(fx, z, 1)) > EPS);

            X = z;
            Y = F.Solve(X);
            E = Math.Abs(Differentiate.Derivative(fx, z, 1));


            a = A;
            b = B;
            int i = 0;

            while (i < It)
            {
                z = (a + b) / 2;

                if (Differentiate.Derivative(fx, z, 1) < 0)
                {
                    a = z;
                }
                else
                {
                    b = z;
                }

                i++;
            } 

            Acc = Math.Abs(Differentiate.Derivative(fx, z, 1));
        }
    }

    internal class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("\n\tSolution of equations with one variable.\n");
            Console.WriteLine("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
            Console.WriteLine("\tEquations:");
            Console.WriteLine("\t1) (x+1)^4 - 2 * x^2\n");
            Console.WriteLine("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

            Console.Write("\tEnter equation number: ");
            int n = Convert.ToInt32(Console.ReadLine());

            Console.Write("\tRange start: ");
            double a = Convert.ToDouble(Console.ReadLine());

            Console.Write("\tRange end: ");
            double b = Convert.ToDouble(Console.ReadLine());

            Console.Write("\tCalculation accuracy: ");
            double eps = Convert.ToDouble(Console.ReadLine());

            Console.Write("\tNumber of iterations: ");
            int it = Convert.ToInt32(Console.ReadLine());

            GoldenSection m1 = new GoldenSection(n, a, b, eps, it);
            Console.WriteLine(m1);

            Dichotomy m2 = new Dichotomy(n, a, b, eps, it);
            Console.WriteLine(m2);

            Newton_Raphson m3 = new Newton_Raphson(n, a, b, eps, it);
            Console.WriteLine(m3);

            MiddlePoints m4 = new MiddlePoints(n, a, b, eps, it);
            Console.WriteLine(m4);


            Console.WriteLine($"\n\n\tx0 - x: \n" +
                $"\tMethod {m1.GetType()} - {m1.F.X - m1.X:f10}\n" +
                $"\tMethod {m2.GetType()} - {m2.F.X - m2.X:f10}\n" +
                $"\tMethod {m3.GetType()} - {m3.F.X - m3.X:f10}\n" +
                $"\tMethod {m4.GetType()} - {m4.F.X - m4.X:f10}\n");

            Console.WriteLine($"\tEps: \n" +
                $"\tMethod {m1.GetType()} - {m1.E:f10}\n" +
                $"\tMethod {m2.GetType()} - {m2.E:f10}\n" +
                $"\tMethod {m3.GetType()} - {m3.E:f10}\n" +
                $"\tMethod {m4.GetType()} - {m4.E:f10}\n");

            Console.WriteLine($"\tNumber of iterations: \n" +
                $"\tMethod {m1.GetType()} - {m1.I:f10}\n" +
                $"\tMethod {m2.GetType()} - {m2.I:f10}\n" +
                $"\tMethod {m3.GetType()} - {m3.I:f10}\n" +
                $"\tMethod {m4.GetType()} - {m4.I:f10}\n");

            Console.WriteLine($"\tAccuracy at {it} iterations: \n" +
                $"\tMethod {m1.GetType()} - {m1.Acc:f10}\n" +
                $"\tMethod {m2.GetType()} - {m2.Acc:f10}\n" +
                $"\tMethod {m3.GetType()} - {m3.Acc:f10}\n" +
                $"\tMethod {m4.GetType()} - {m4.Acc:f10}\n");
        }
    }
}
