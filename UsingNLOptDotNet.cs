// Implements sample problem from http://ab-initio.mit.edu/wiki/index.php/NLopt_Tutorial. See this wiki for
// details of Prof. Steven Johnson's nlopt. 

using System;
using System.Collections.Generic;
using System.Text;

// Uses NLOptDotNet wrapper, written in C++/CLI.
 
// Add a reference to NLOptDotNet.dll (set it's 'copy local' property to 'True' so that a copy of the dll is present in your 
// application's current directory), and make sure libnlopt-0.dll is also in your application's current directory at execution. 
// For resons that I'm not 100% clear on, but are can likely be relaxed with 'include' directories or similar, both these dll's 
// must be in the application's execution directory at runtime.

using NLOptDotNet; // This namespace contains a managed wrapper to the unamanaged libnlopt-0.dll

namespace UsingNLOptDotNet
{
    class Program
    {
        // verbose console output trigger
        static bool VERBOSE_OUTPUT = false;
        
        [STAThread]
        static void Main(string[] args)
        {
            // The optimization algorithm
            Algorithm main_alg = Algorithm.AUGLAG;
            // The local/subsidiary optimization algorithm, which may or may not be used (see NLOpt wiki)
            Algorithm secondary_alg = Algorithm.GN_CRS2_LM;

            // Create optimization object, setting algorithm type and number of variables
            NLOptWrapper wrapper = new NLOptWrapper(main_alg, 2);
            
            // Turn verbose (console) output on or off
            wrapper.VERBOSE_OUTPUT = VERBOSE_OUTPUT;

            // Set some stopping criteria
            wrapper.MaxTime = 100; //in seconds
            wrapper.XTolRelative = 1e-4;
            wrapper.FTolAbsolute = 1e-2;
            wrapper.MaxEval = 5000;
            wrapper.StopVal = 0.55;
                        
            // Create design vector and set initial guess 
            double[] x = new double[2] { 1.234, 5.678 };
            
            // Apply lower bounds on the design vector
            wrapper.SetUpperBounds(new double[2] { 10000, 10000 });
            wrapper.SetLowerBounds(new double[2] { -10000, 0.0001 });
                        
            // Create local optimization object, if appropriate (all AUGLAG formulations and MLSL require this)
            // ...Note that this has to be done before the main wrapper's objectives or constraints are set
            // to prevent crashing. Don't know exactly why yet. Part of the joys of mixing unmanaged and managed code!               
            if (main_alg.ToString().Contains("AUGLAG") || main_alg.ToString().Contains("MLSL"))
            {
                NLOptWrapper local_wrapper = new NLOptWrapper(secondary_alg, wrapper.Dimension);
                /* add stopping criteria */ local_wrapper.XTolRelative = 1e-4; local_wrapper.FTolRelative = 1e-2;
                
                // Haven't figured out whether the local or main algorithm stopping criteria dominate, and when. 
                // Need to inspect nlopt source code to be sure of specifics, and differences between AUGLAG and MLSL
                local_wrapper.VERBOSE_OUTPUT = VERBOSE_OUTPUT;
                wrapper.SetLocalOptimizer(local_wrapper);
            }

            // Delegate the objective function. Data can be passed in as type System.Object and cast in the objective function
            wrapper.SetMinObjective(new FunctionDelegate(Objective), null);
            
            // Add inequality constraints. Data can be passed in as type System.Object. Feasibility tolerance passed as an argument
            my_constraint_data[] data = new my_constraint_data[2] { new my_constraint_data(2, 0), new my_constraint_data(-1, 1) };
            wrapper.AddInequalityConstraint(new FunctionDelegate(Constraint), data[0], 1e-8);
            wrapper.AddInequalityConstraint(new FunctionDelegate(Constraint), data[1], 1e-8);
                        
            // create variable to store min objective
            double minf = 0;

            //Run the optimization, passing minf by reference. NLOptDotNet.Result returned
            Result r = wrapper.Optimize(x, ref minf);

            Console.WriteLine("\nFound minimum after " + neval + " objective function evaluations");
            Console.WriteLine("Found minimum at f(" + x[0] + ", " + x[1] + ") = " + minf);
            Console.WriteLine("\nNLOpt Result: " + r.ToString());
        }
                
        // Function signature for objective and constraint function enforced by the delegate 'NLOptDotNet.FunctionDelegate'
        static int neval = 0; // counts how many times the objective fnction is called        
        static double Objective(double[] x, ref double[] grad, Object data)
        {
            Console.WriteLine(neval++);
            if (VERBOSE_OUTPUT) { Console.WriteLine(".NET Objective #" + neval + ": objective function called with point (" + x[0] + ", " + x[1] + ")"); }
            if (grad != null)
            {// NLOptDotNet.NLOptWrapper will send a gradient vector by ref when the algorithm used is gradient-based
                grad[0] = 0.0;
                grad[1] = 0.5 / Math.Sqrt(x[1]);
                if (VERBOSE_OUTPUT) { Console.WriteLine(".NET Objective: Gradient on objective function is (" + grad[0] + ", " + grad[1] + ")"); }
            }
            return Math.Sqrt(x[1]);
        }
        static double Constraint(double[] x, ref double[] grad, Object data)
        {
            if (VERBOSE_OUTPUT) { Console.WriteLine(".NET Constraint: constraint function called with point (" + x[0] + ", " + x[1] + ")"); }
            // Data passed as System.Object to parameterize the constraint functions. Convert to my_constraint_data struct
            my_constraint_data d = (my_constraint_data)data;
            double a = d.a, b = d.b;
            if (grad != null)
            {
                grad[0] = 3 * a * (a * x[0] + b) * (a * x[0] + b);
                grad[1] = -1.0;
                if (VERBOSE_OUTPUT) { Console.WriteLine(".NET Constraint: Gradient on constraint function is (" + grad[0] + ", " + grad[1] + ")"); }
            }
            double answer = ((a * x[0] + b) * (a * x[0] + b) * (a * x[0] + b) - x[1]);
            return answer;
        }
        // A struct to pass data to constraint functions
        public struct my_constraint_data
        {
            public double a;
            public double b;
            public my_constraint_data(double a, double b)
            {
                this.a = a; this.b = b;
            }
        }
    }
}
