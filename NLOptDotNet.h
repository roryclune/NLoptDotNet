// NLOptDotNet.h
//		Functionality provided by nlopt.hpp that I HAVE NOT yet wrapped:
//		vector-valued constraints, bound constraints passed by reference, version number

#pragma once
using namespace System; //CLR
using namespace System::Collections::Generic;
using namespace System::Runtime::InteropServices;

using namespace std;	//old-school
#include <algorithm>	//for std::copy
#include <exception>

#include "nlopt.hpp"	//to be wrapped
using namespace nlopt;

namespace NLOptDotNet {
	
	//-------------------------------------------------------------
	// The types of algorithm that can be used.
	// Need to cast betwen this and native unmananged enum whenever used
	//-------------------------------------------------------------	
	public enum class Algorithm
	{
		GN_DIRECT = 0,
		GN_DIRECT_L, GN_DIRECT_L_RAND, GN_DIRECT_NOSCAL, GN_DIRECT_L_NOSCAL, GN_DIRECT_L_RAND_NOSCAL, GN_ORIG_DIRECT,
		GN_ORIG_DIRECT_L, GD_STOGO, GD_STOGO_RAND, LD_LBFGS_NOCEDAL, LD_LBFGS, LN_PRAXIS, LD_VAR1, LD_VAR2, LD_TNEWTON,
		LD_TNEWTON_RESTART, LD_TNEWTON_PRECOND, LD_TNEWTON_PRECOND_RESTART, GN_CRS2_LM, GN_MLSL, GD_MLSL, GN_MLSL_LDS,
		GD_MLSL_LDS, LD_MMA, LN_COBYLA, LN_NEWUOA, LN_NEWUOA_BOUND, LN_NELDERMEAD, LN_SBPLX, LN_AUGLAG, LD_AUGLAG,
		LN_AUGLAG_EQ, LD_AUGLAG_EQ, LN_BOBYQA, GN_ISRES, AUGLAG, AUGLAG_EQ, G_MLSL, G_MLSL_LDS, LD_SLSQP,
		NUM_ALGORITHMS /* not an algorithm, just the number of them */
	};
	
	//-------------------------------------------------------------
	// Possible nlopt results. 
	// Need to cast betwen this and native unmananged enum whenever used
	//-------------------------------------------------------------	
	public enum class Result {
		FAILURE = -1, /* generic failure code */
		INVALID_ARGS = -2, OUT_OF_MEMORY = -3, ROUNDOFF_LIMITED = -4, FORCED_STOP = -5, SUCCESS = 1, /* generic success code */
		STOPVAL_REACHED = 2, FTOL_REACHED = 3, XTOL_REACHED = 4, MAXEVAL_REACHED = 5, MAXTIME_REACHED = 6 
	};
	
	//-------------------------------------------------------------
	// The different types of function that can be delegated and used
	// Allows the use of one private 'vfuncWrapper' method for multiple funciton types
	//-------------------------------------------------------------	
	private enum class FunctionType
	{
		MINIMUM_OBJECTIVE = 1, MAXIMUM_OBJECTIVE = 2, INEQUALITY_CONSTRAINT = 3, EQUALITY_CONSTRAINT = 4
	};
	
	//-------------------------------------------------------------
	// A delegate, used by .NET in place of function pointers
	// Its parameters have to be manually converted to native types before being passed to nlopt
	//-------------------------------------------------------------
	public delegate double FunctionDelegate(array<double>^ x, array<double>^% grad, Object^ data);
	
	//-------------------------------------------------------------	
	// The function_data struct wraps whatever data object the client supplies, and adds an id tag
	// so the wrapper functions know which of the .NET delegates it should call when processing
	// a native callback
	//-------------------------------------------------------------	
	private ref struct function_data
	{ 
		public: 
			int id; 
			Object^ data;
	};

	//-------------------------------------------------------------
	// The main wrapper class - OptWrapper. Ideally this should contain only a pointer to the native class and methods 
	// that wrap native methods. I have had to include some data members to manage conversions from managed to unamanged
	// types and memory/Garbage Collector issues
	//-------------------------------------------------------------
	public ref class NLOptWrapper
	{	
	private:
		// a lock token to explore thread safety
		static bool locked = false;

		//a C++ pointer to an nlopt::opt object
		nlopt::opt* _opt;
		// instances of the delegates for the .NET functions
		List<FunctionDelegate^>^ _funcDelegatesList;
		// A counter that tracks how many delegates (constraints or objectives) have been passed in by the client
		int nfunc;
		// verbose output setting
		bool _verbose_output;

		// GCHandles are used to prevent delegates passed (by void pointers)to unmanaged code from being relocated on 
		// the managed heap by the .NET garbage collector
		List<GCHandle>^ _GCHandles;

		// A private delegate to handle the callback from libnlopt-0.dll. Will convert the parameters and call the
		// .NET function, and then convert back any parameters that were expected to be passed by reference
		delegate double vfuncWrapperDel(const std::vector<double> &x, std::vector<double> &grad, void *data);
		
		//-------------------------------------------------------------		
		double vfuncWrapper(const std::vector<double> &x, std::vector<double> &grad, void *data)
			//-------------------------------------------------------------		
		{
			if(_verbose_output) { Console::WriteLine("\nin NLOptWrapper.vfuncWrapper(): callback to wrapper made by libnlopt-0.dll"); }
			// Copy the native variables into the .NET variables
			int n = x.size();			
			array<double>^ mx = gcnew array<double>(n);						
			for(int i = 0; i<n ; i++)
			{
				mx[i] = x[i]; //can't use Marshal::Copy here because of the const modifier on parameter x			
			}			
			// Copy the native grad vector into managed .NET, if the vector exists
			array<double>^ mGrad;
			if(!grad.empty())
			{
				mGrad = gcnew array<double>(n);
				Marshal::Copy(IntPtr(&grad), mGrad, 0, n);				
			}

			//Recover the function_data struct from the native void pointer, which points to the managed heap
			GCHandle h = GCHandle::FromIntPtr(IntPtr(data));
			function_data^ mData = (function_data^)h.Target;			

			// Call the .NET method through its delegate, which is in the id'th location in the relevant List
			FunctionDelegate^ fd = _funcDelegatesList[mData->id];
			double answer = fd(mx, mGrad, mData->data);

			if(_verbose_output) { Console::WriteLine("in NLOptWrapper.vfuncWrapper(): Value returned to nlopt = " + answer);}
			if(!grad.empty())
			{
				//Copy the variables that were expected to have been passed by ref by libnlopt-0.dll (grad vector only)
				Marshal::Copy(mGrad, 0, IntPtr(&grad[0]), mGrad->Length);			
			}
			return answer;
		}
		void AddFunction(FunctionType ftype, FunctionDelegate^ funcDelegate, Object^ data, double tol)
		{
			//Store the delegate just passed in so that it remains in scope after this method terminates 
			_funcDelegatesList->Add(funcDelegate);

			// Create a function pointer from a delegate to vfuncWrapper. Note that we forget all about the .NET delegate 
			// for now. It just gets added to a list and retrieved based on its function_data struct in the vfuncWrapper method
			vfuncWrapperDel^ fp = gcnew vfuncWrapperDel(this, &NLOptWrapper::vfuncWrapper);
			_GCHandles->Add(GCHandle::Alloc(fp));
			IntPtr ip = Marshal::GetFunctionPointerForDelegate(fp);
			nlopt::vfunc cb = static_cast<nlopt::vfunc>(ip.ToPointer());

			// create a function_data struct, wrapping the data passed in and giving it an identifier to retrieve the 
			// .NET delegate later
			function_data^ fdata = gcnew function_data();
			fdata->id = nfunc; fdata->data = data;
			nfunc++;
			// A GCHandle on fdata, the function_data struct, maintains its location on the managed heap after it leaves scope
			GCHandle gch = GCHandle::Alloc(fdata);
			_GCHandles->Add(gch);
			void* voidPtr = GCHandle::ToIntPtr(gch).ToPointer();

			switch(ftype)
			{
			case FunctionType::MINIMUM_OBJECTIVE:
				_opt->set_min_objective(cb, voidPtr);
				break;
			case FunctionType::MAXIMUM_OBJECTIVE:
				_opt->set_max_objective(cb, voidPtr);
				break;
			case FunctionType::EQUALITY_CONSTRAINT:
				_opt->add_equality_constraint(cb, voidPtr, tol);
				break;			
			case FunctionType::INEQUALITY_CONSTRAINT:				
				_opt->add_inequality_constraint(cb, voidPtr, tol);
				break;
			default:break;						
			}
		}
		//-------------------------------------------------------------
		property nlopt::opt* NativeOptPointer
			//-------------------------------------------------------------		
		{
			//Exposes the _opt object
			nlopt::opt* get()
			{ 
				return _opt;			
			}
		}

	public:
		//-------------------------------------------------------------
		NLOptWrapper(){} // The mandatory default constructor			
		NLOptWrapper(Algorithm alg, int nvar)
		//-------------------------------------------------------------
		{
			nlopt::algorithm unmanaged_alg = static_cast<nlopt::algorithm>(alg);			
			_opt = new nlopt::opt(unmanaged_alg, nvar);
			_funcDelegatesList = gcnew List<FunctionDelegate^>();
			_GCHandles = gcnew List<GCHandle>();
			nfunc = 0;
			_verbose_output = false;
		}
		//-------------------------------------------------------------
		~NLOptWrapper()
			//-------------------------------------------------------------
		{
			try
			{
				// cleanup managed resources
				for each (GCHandle gch in _GCHandles) { gch.Free(); }
				delete _GCHandles; _GCHandles = nullptr;				
				delete _funcDelegatesList; _funcDelegatesList = nullptr;			
			}
			catch(System::Exception^ e)
			{
				Console::WriteLine("Cleanup of managed resources failed for some reason");			
			}
			finally
			{
				// cleanup native resources even if member variables could not be disposed
				nfunc = 0;
				delete _opt; _opt = 0; // Delete the pointer to the nlopt::opt objective				
			}
			// Free all the GCHandle objects that are pinning function delegates and structs, avoiding memory leak
			if(_verbose_output){Console::WriteLine("NLOptWrapper's destructor run. All GCHandles freed and member variables deleted. Managed resources cleaned up");}
		}
		//-------------------------------------------------------------
		void SetMinObjective(FunctionDelegate^ objFuncDelegate, Object^ data)
			//-------------------------------------------------------------		
			//WRAPS...void set_min_objective(vfunc vf, void *f_data)
		{
			AddFunction(FunctionType::MINIMUM_OBJECTIVE, objFuncDelegate, data, 0);			
		}
		void SetMaxObjective(FunctionDelegate^ objFuncDelegate, Object^ data)
			//-------------------------------------------------------------		
			//WRAPS...void set_max_objective(vfunc vf, void *f_data)
		{
			AddFunction(FunctionType::MAXIMUM_OBJECTIVE, objFuncDelegate, data, 0);			
		}
		//-------------------------------------------------------------
		void AddInequalityConstraint(FunctionDelegate^ conFuncDelegate, Object^ data, double tol)
			//-------------------------------------------------------------
			//WRAPS...void add_inequality_constraint(vfunc vf, void *f_data, double tol=0)			
		{
			AddFunction(FunctionType::INEQUALITY_CONSTRAINT, conFuncDelegate, data, tol);
		}
		void AddEqualityConstraint(FunctionDelegate^ conFuncDelegate, Object^ data, double tol)
			//-------------------------------------------------------------
			//WRAPS...void add_inequality_constraint(vfunc vf, void *f_data, double tol=0)			
		{
			AddFunction(FunctionType::EQUALITY_CONSTRAINT, conFuncDelegate, data, tol);
		}
		//------------------------------------------------------------		
		void SetLocalOptimizer(NLOptWrapper^ LocalOpt)
		//	-------------------------------------------------------------		
		{
			// WRAPS...void nlopt::opt::set_local_optimizer(const nlopt::opt &local_opt)		
			// from <http://ab-initio.mit.edu> nlopt::opt object whose parameters are used to determine the local search algorithm, 
			// its stopping criteria, and other algorithm parameters. (However, the objective function, bounds, and nonlinear-
			// constraint parameters of local_opt are ignored.) The dimension n of local_opt must match that of opt.
			
			_GCHandles->Add(GCHandle::Alloc(LocalOpt));			// to preserve location in memory
			nlopt::opt* _sub_opt = LocalOpt->NativeOptPointer;	// returns a pointer to the object
			nlopt::opt test = *_sub_opt;
			_opt->set_local_optimizer(*_sub_opt);			
		}		
		//-------------------------------------------------------------
		Result Optimize(array<double>^x, double% f) //both passed by reference
			//-------------------------------------------------------------					
		{			
			// First, create native data types 
			vector<double> native_x(x->Length);
			Marshal::Copy(x, 0, IntPtr(&native_x[0]), x->Length);
			double native_f = 0; //will be passed by reference

			// Run the optimization, passing the new native variables by reference
			nlopt::result r;
			try
			{
				if(_verbose_output){Console::WriteLine("in NLOptWrapper.Optimize(): starting point about to be sent to nlopt dll by wrapper");}
				r = _opt->optimize(native_x, native_f);
			}
			// Catch exceptions
			// Below are the specific types of c++ (std::) exception that can be thrown
			catch(bad_alloc& e)
			{
				Console::WriteLine("\n------\n\"bad_alloc\" thrown during NLopt's Optimize() method");
				Console::WriteLine("Ran out of memory (a memory allocation failed). \n'what()' message = " + e.what()->ToString());
				Console::WriteLine("------\n");
				r = OUT_OF_MEMORY;
			}
			catch(invalid_argument& e)
			{
				Console::WriteLine("\n------\n\"Invalid arguments\" thrown during NLopt's optimize() method");
				Console::WriteLine("(e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etc.). \n'what()' message = " + e.what()->ToString());				
				Console::WriteLine("------\n");
				r = INVALID_ARGS;
			}
			catch(roundoff_limited& e)
			{
				Console::WriteLine("\n------\n\"Roundoff limited\" thrown during NLopt's optimize() method");
				Console::WriteLine("Halted because roundoff errors limited progress. \n'what()' message = " + e.what()->ToString());
				Console::WriteLine("------\n");
				r = ROUNDOFF_LIMITED;
			}
			catch(forced_stop& e)
			{
				Console::WriteLine("\n------\n\"Forced stop\" exception thrown during NLopt's optimize() method");
				Console::WriteLine("Halted because of a forced termination: the user called nlopt::opt::force_stop() from the user’s objective function or threw an nlopt::forced_stop exception. 'what()' message = " + e.what()->ToString());
				Console::WriteLine("------\n");
				r = FORCED_STOP;
			}
			catch(runtime_error& e)
			{
				Console::WriteLine("\n------\nNLOpt exception thrown during NLopt's optimize() method");
				Console::WriteLine("Generic failure. 'what()' message = " + e.what()->ToString());
				Console::WriteLine("------\n");
				r = FAILURE;			
			}
			// And finally, a generic catch:
			catch(System::Exception^ e)
			{
				Console::WriteLine("\n-------\nNLOpt threw an exception during optimize() method\n");
				while(e)
				{
					Console::WriteLine("\nSource: " + e->Source);
					Console::WriteLine("\nMessage: " + e->Message);
					Console::WriteLine("\nStack Trace: " + e->StackTrace);			
					Console::WriteLine("\n-------\n");
					e = e->InnerException;
					r = FAILURE;
				}
			}			
			if(_verbose_output){Console::WriteLine("OptWrapper.Optimize(): after optimization");}

			// Copy the native variables into the (by reference) .NET variables
			Marshal::Copy(IntPtr(&native_x[0]), x, 0, native_x.size());
			f = native_f;
			// Case the result from native to managed enum and return
			return static_cast<Result>(r);			
		}
		//-------------------------------------------------------------				
		void Force_Stop()
		{
			_opt->force_stop();			
		}
		//-------------------------------------------------------------				
		property Algorithm GetAlgorithm
			//-------------------------------------------------------------				
		{
			Algorithm get(){ return (Algorithm)_opt->get_algorithm(); }			
		}
		//-------------------------------------------------------------				
		property int Dimension
			//-------------------------------------------------------------
		{
			int get(){ return _opt->get_dimension(); }
		}
		//-------------------------------------------------------------				
		property bool VERBOSE_OUTPUT
			//-------------------------------------------------------------				
		{
			bool get(){ return _verbose_output; }
			void set(bool value){_verbose_output = value;}
		}		
		//-------------------------------------------------------------				
		property double StopVal
			//-------------------------------------------------------------				
		{
			// Stop when an objectve value of at least stopval is reached
			double get(){ return _opt->get_stopval(); }
			void set(double stopval){_opt->set_stopval(stopval);}			
		}
		//-------------------------------------------------------------				
		property double FTolRelative
			//-------------------------------------------------------------				
		{
			// Set relative tolerance on function value.
			double get(){ return _opt->get_ftol_rel(); }
			void set(double tol){_opt->set_ftol_rel(tol);}
		}
		//-------------------------------------------------------------				
		property double FTolAbsolute
			//-------------------------------------------------------------				
		{
			// Set absolute tolerance on function value.
			double get(){ return _opt->get_ftol_abs(); }
			void set(double tol){_opt->set_ftol_abs(tol);}			
		}
		//-------------------------------------------------------------				
		property double XTolRelative
			//-------------------------------------------------------------				
		{
			// Set relative tolerance on optimization parameters (design vector)
			double get(){ return _opt->get_xtol_rel(); }
			void set(double tol){_opt->set_xtol_rel(tol);}
		}
		//-------------------------------------------------------------				
		property array<double>^ XTolAbsolute
			//-------------------------------------------------------------				
		{
			// Set absolute tolerance on optimization parameters (design vector)
			array<double>^ get()
			{
				array<double>^ tols = gcnew array<double>(Dimension);				
				std::vector<double> native(tols->Length);
				_opt->get_xtol_abs(native);
				Marshal::Copy(IntPtr(&native), tols, 0, native.size());
				return tols;			
			}
			void set(array<double>^ tols)
			{
				std::vector<double> native(tols->Length);
				Marshal::Copy(tols, 0, IntPtr(&native[0]), tols->Length);
				_opt->set_xtol_abs(native);
			}
		}
		//-------------------------------------------------------------				
		property int MaxEval
			//-------------------------------------------------------------				
		{
			// Stop when the number of function evaluations exceeds maxeval
			int get(){ return _opt->get_maxeval(); }
			void set(int maxeval){_opt->set_maxeval(maxeval);}
		}
		//-------------------------------------------------------------				
		property double MaxTime
			//-------------------------------------------------------------				
		{
			// Stop when the optimization time (in seconds) exceeds maxtime
			double get(){ return _opt->get_maxtime(); }
			void set(double seconds){_opt->set_maxtime(seconds);}
		}
		//-------------------------------------------------------------				
		property array<double>^ InitialStepSize
			//-------------------------------------------------------------				
		{
			// For derivative-free algorithms. Size of initial perturbation
			array<double>^ get()
			{ 
				std::vector<double> steps = _opt->get_lower_bounds();//CHANGE//
				array<double>^ dx = gcnew array<double>(steps.size());			
				Marshal::Copy(IntPtr(&steps), dx, 0, steps.size());
				return dx;
			}
			void set(array<double>^ dx)
			{
				std::vector<double> steps(dx->Length);
				Marshal::Copy(dx, 0, IntPtr(&steps[0]), dx->Length);
				_opt->set_initial_step(steps);
			}
		}
		//-------------------------------------------------------------				
		array<double>^ GetInitialStepSize(array<double>^ initial_x)
			//-------------------------------------------------------------				
		{
			// For derivative-free algorithms. Size of initial perturbation. Starting point has to be passed 
			// in since perturbation depends on this if client has not already set it

			std::vector<double> native_x(initial_x->Length);
			Marshal::Copy(initial_x, 0, IntPtr(&native_x[0]), initial_x->Length);

			std::vector<double> native_dx(initial_x->Length);			
			_opt->get_initial_step(native_x, native_dx);

			array<double>^ dx = gcnew array<double>(Dimension);			
			Marshal::Copy(IntPtr(&native_dx), dx, 0, native_dx.size());
			return dx;
		}
		//-------------------------------------------------------------						
		void SetInitialStepSize(array<double>^ dx)
			//-------------------------------------------------------------		
		{
			std::vector<double> native(dx->Length);
			Marshal::Copy(dx, 0, IntPtr(&native[0]), dx->Length);
			_opt->set_initial_step(native);
		}
		//-------------------------------------------------------------
		void SetLowerBounds(array<double>^ lb)
			//-------------------------------------------------------------
		{	
			std::vector<double> lower(lb->Length);
			Marshal::Copy(lb, 0, IntPtr(&lower[0]), lb->Length);			
			_opt->set_lower_bounds(lower);			
		}
		//-------------------------------------------------------------
		void SetLowerBounds(double lb)
			//-------------------------------------------------------------
		{
			_opt->set_lower_bounds(std::vector<double>(Dimension, lb));						
		}
		//-------------------------------------------------------------
		array<double>^ GetLowerBounds()
			//-------------------------------------------------------------
		{
			std::vector<double> lower = _opt->get_lower_bounds();
			array<double>^ lb = gcnew array<double>(lower.size());			
			Marshal::Copy(IntPtr(&lower), lb, 0, lower.size());
			return lb;
		}
		//-------------------------------------------------------------
		void SetUpperBounds(array<double>^ ub)
			//-------------------------------------------------------------
		{
			std::vector<double> upper(ub->Length);
			Marshal::Copy(ub, 0, IntPtr(&upper[0]), ub->Length);
			_opt->set_upper_bounds(upper);			
		}
		//-------------------------------------------------------------
		void SetUpperBounds(double ub)
			//-------------------------------------------------------------
		{
			_opt->set_upper_bounds(std::vector<double>(Dimension, ub));			
		}
		//-------------------------------------------------------------
		array<double>^ GetUpperBounds()
			//-------------------------------------------------------------
		{
			std::vector<double> upper = _opt->get_upper_bounds();
			array<double>^ ub = gcnew array<double>(upper.size());
			Marshal::Copy(IntPtr(&upper), ub, 0, upper.size());			
			return ub;
		}
		//-------------------------------------------------------------		
		void SetPseudoRandomSeed(int seed)
			//-------------------------------------------------------------		
		{
			//For stochastic algorithms, instead ofusing system time as seed
			nlopt::srand(seed);
		}
		//-------------------------------------------------------------		
		void PseudoRandomSeed_Time(int seed)
			//-------------------------------------------------------------		
		{
			//For stochastic algorithms, reset seed to system time
			nlopt::srand_time();
		}
		//-------------------------------------------------------------
		property int VectorStorage
			//-------------------------------------------------------------				
		{
			//For limited memory quasi-newton algorithms. Number M of stored vectors
			int get(){return _opt->get_vector_storage();}
			void set(int value)
			{
				_opt->set_vector_storage(value);				
			}		
		}		
	};
}
