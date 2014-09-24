#ifndef _ARRAY_FUNCTIONS_H_
#define _ARRAY_FUNCTIONS_H_

// Sorta provides intermediate wrapper functions for many of the
// operations; it's just easier to do this with templating here than
// in cython.

#include <ilconcert/iloalg.h>
#include <ilconcert/iloenv.h>
#include <ilconcert/ilosolution.h>
#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplexi.h>
#include <sstream>

#include "debug.h"
#include "optimizations.h"
#include "simple_shared_ptr.h"
#include "containers.hpp"
#include "operators.hpp"

using namespace std;

#define MODEL_UNBOUNDED 2
#define MODEL_INFEASABLE 3
#define MODEL_UNBOUNDED_OR_INFEASABLE 4

#define CPX_STAT_OPTIMAL 1
#define CPX_STAT_UNBOUNDED 2
#define CPX_STAT_INFEASIBLE 3
#define CPX_STAT_INForUNBD 4
#define CPX_STAT_OPTIMAL_INFEAS 5
#define CPX_STAT_NUM_BEST 6
#define CPX_STAT_ABORT_IT_LIM 10
#define CPX_STAT_ABORT_TIME_LIM 11
#define CPX_STAT_ABORT_OBJ_LIM 12
#define CPX_STAT_ABORT_USER 13
#define CPX_STAT_OPTIMAL_FACE_UNBOUNDED 20
#define CPX_STAT_ABORT_PRIM_OBJ_LIM 21
#define CPX_STAT_ABORT_DUAL_OBJ_LIM 22
#define CPXMIP_OPTIMAL 101
#define CPXMIP_OPTIMAL_TOL 102
#define CPXMIP_INFEASIBLE 103
#define CPXMIP_SOL_LIM 104
#define CPXMIP_NODE_LIM_FEAS 105
#define CPXMIP_NODE_LIM_INFEAS 106
#define CPXMIP_TIME_LIM_FEAS 107
#define CPXMIP_TIME_LIM_INFEAS 108
#define CPXMIP_FAIL_FEAS 109
#define CPXMIP_FAIL_INFEAS 110
#define CPXMIP_MEM_LIM_FEAS 111
#define CPXMIP_MEM_LIM_INFEAS 112
#define CPXMIP_ABORT_FEAS 113
#define CPXMIP_ABORT_INFEAS 114
#define CPXMIP_OPTIMAL_INFEAS 115
#define CPXMIP_FAIL_FEAS_NO_TREE 116
#define CPXMIP_FAIL_INFEAS_NO_TREE 117
#define CPXMIP_UNBOUNDED 118
#define CPXMIP_INForUNBD 119

////////////////////////////////////////////////////////////////////////////////
// Now the same for the model

class CPlexModelInterface {
public:
  struct Status {
    Status(const char* _message, int _error_code = 1)
      : error_code(_error_code), message(_message)
    {
    }

    Status()
      : error_code(0), message("Okay.")
    {
    }
		
    int error_code;
    const char* message;
  };
  
  CPlexModelInterface(IloEnv _env) 
    : env(_env), model(_env), solver(_env), current_objective(NULL), model_extracted(false), model_solved(false)
  {
  }

  Status addVariables(const ExpressionArray& expr)
  {
    try{
      model.add(expr.variables());
    } catch(IloException& e) {
      return Status(e.getMessage());
    }

    return Status();
  }

  Status addConstraint(const ConstraintArray& cstr)
  {
    model_solved = false;

    try{
      model.add(cstr.constraint());
    } catch(IloException& e) {
      return Status(e.getMessage());
    }

    return Status();
  }

  Status setObjective(const ExpressionArray& expr, bool maximize)
  {
    model_solved = false;
	    
    try {
      if(current_objective != NULL) {
	model.remove(*current_objective);
	delete current_objective;
	current_objective = NULL;
      }
    } catch(IloException& e) {
      return Status(e.getMessage());
    }

    if(expr.shape(0) != 1 || expr.shape(1) != 1)
      return Status("Objective must be a scalar (or 1x1 matrix) expression.");

    try {
		
      if(maximize)
	current_objective = new IloObjective(IloMaximize(env, expr(0,0)));
      else
	current_objective = new IloObjective(IloMinimize(env, expr(0,0)));
		
      model.add(*current_objective);
    }
    catch(IloException& e) {
      return Status(e.getMessage());
    }

    return Status();
  }
    
  Status removeConstraint(const ConstraintArray& csr)
  {
    model_solved = false;

    try {
      model.remove(csr.constraint());
    } catch(IloException& e) {
      return Status(e.getMessage());
    }

    return Status();
  }

  Status setStartingValues(const ExpressionArray& expr, const NumericalArray& numr)
  {
    try {
      if(!model_extracted)
	extractModel();

      if(!expr.hasVar())
	return Status("Only variables may be set, not expressions.");

      IloNumArray X(env, expr.size());

      if(expr.isComplete()) 
	{
	  for(long i = 0; i < expr.shape(0); ++i)
	    for(long j = 0; j < expr.shape(1); ++j)
	      X[expr.getIndex(i,j)] = numr(i,j);

	  solver.getImpl()->setVectors(X, 0, expr.variables(), 0, 0, 0);
	} 
      else
	{
	  IloNumVarArray xpr(env, expr.size());

	  for(long i = 0; i < expr.shape(0); ++i) {
	    for(long j = 0; j < expr.shape(1); ++j)
	      {
		const IloNumVarArray& vars = expr.variables();
			    
		X[i*expr.shape(1) + j] = numr(i,j);
		xpr[i*expr.shape(1) + j] = vars[expr.getIndex(i,j)];
	      }
	  }
	  solver.getImpl()->setVectors(X, 0, xpr, 0, 0, 0);
	}
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }

    return Status();
  }
    
  template <typename Param, typename V>
  Status setParameter(const Param& p, V value)
  {
    try{
      solver.setParam(p, value);
    } catch(IloException& e){
      return Status(e.getMessage());
    }

    return Status();
  }

  Status handleCplexStatus(IloCplex::CplexStatus status) 
  {
    switch(status) {
      /*
    case IloCplex::Optimal:
      return Status("Optimal solution is available", CPX_STAT_OPTIMAL);
    case IloCplex::Unbounded:
      return Status("Model has an Unbounded ray", CPX_STAT_UNBOUNDED);
    case IloCplex::Infeasible:
      return Status("Model is proved Infeasible", CPX_STAT_INFEASIBLE);
    case IloCplex::InfOrUnbd:
      return Status("Model is proved either Infeasible or Unbounded", CPX_STAT_INForUNBD);
    case IloCplex::OptimalInfeas:
      return Status("Optimal solution is available, but with infeasibilities after unscaling", CPX_STAT_OPTIMAL_INFEAS);
      */
    case IloCplex::NumBest:
      return Status("Solution is available, but not proved optimal, due to numerical difficulties during optimization", CPX_STAT_NUM_BEST);
    case IloCplex::AbortItLim:
      return Status("Aborted due to an iteration limit", CPX_STAT_ABORT_IT_LIM);
      /*
    case IloCplex::AbortTimeLim:
      return Status("Aborted due to a time limit", CPX_STAT_ABORT_TIME_LIM);
      */
    case IloCplex::AbortObjLim:
      return Status("Aborted due to an objective limit", CPX_STAT_ABORT_OBJ_LIM);
    case IloCplex::AbortUser:
      return Status("Aborted on user request", CPX_STAT_ABORT_USER);
    case IloCplex::OptimalFaceUnbounded:
      return Status("Model has Unbounded optimal face", CPX_STAT_OPTIMAL_FACE_UNBOUNDED);
    case IloCplex::AbortPrimObjLim:
      return Status("Aborted due to a primal obj limit", CPX_STAT_ABORT_PRIM_OBJ_LIM);
    case IloCplex::AbortDualObjLim:
      return Status("Aborted due to a dual obj limit", CPX_STAT_ABORT_DUAL_OBJ_LIM);  
    case IloCplex::Optimal:
      return Status();
      //"Optimal integer solution found", CPXMIP_OPTIMAL);
    case IloCplex::OptimalTol:
      return Status();
      //"Optimal sol. within epgap or epagap tolerance found", CPXMIP_OPTIMAL_TOL);
    case IloCplex::Infeasible:
      return Status("Integer infeasible", CPXMIP_INFEASIBLE);
    case IloCplex::SolLim:
      return Status("Mixed integer solutions limit exceeded", CPXMIP_SOL_LIM);
    case IloCplex::NodeLimFeas:
      return Status("Node limit exceeded, integer solution exists", CPXMIP_NODE_LIM_FEAS);
    case IloCplex::NodeLimInfeas:
      return Status("Node limit exceeded, no integer solution", CPXMIP_NODE_LIM_INFEAS);
    case IloCplex::AbortTimeLim:
      return Status("Time limit exceeded, integer solution exists", CPXMIP_TIME_LIM_FEAS);
      //case IloCplex::TimeLimInfeas:
      //return Status("Time limit exceeded, no integer solution", CPXMIP_TIME_LIM_INFEAS);
    case IloCplex::FailFeas:
      return Status("Error termination, integer solution exists", CPXMIP_FAIL_FEAS);
    case IloCplex::FailInfeas:
      return Status("Error termination, no integer solution", CPXMIP_FAIL_INFEAS);
    case IloCplex::MemLimFeas:
      return Status("Treememory limit, integer solution exists", CPXMIP_MEM_LIM_FEAS);
    case IloCplex::MemLimInfeas:
      return Status("Treememory limit, no integer solution exists", CPXMIP_MEM_LIM_INFEAS);
      //case IloCplex::AbortFeas:
      //return Status("Aborted, integer solution exists", CPXMIP_ABORT_FEAS);
      //case IloCplex::AbortInfeas:
      //return Status("Aborted, no integer solution", CPXMIP_ABORT_INFEAS);
    case IloCplex::OptimalInfeas:
      return Status("Problem optimal with unscaled infeasibilities", CPXMIP_OPTIMAL_INFEAS);
    case IloCplex::FailFeasNoTree:
      return Status("Out of memory, no tree, integer solution exists", CPXMIP_FAIL_FEAS_NO_TREE);
    case IloCplex::FailInfeasNoTree:
      return Status("Out of memory, no tree, no integer solution", CPXMIP_FAIL_INFEAS_NO_TREE);
    case IloCplex::Unbounded:
      return Status("Model has an Unbounded ray", CPXMIP_UNBOUNDED);
    case IloCplex::InfOrUnbd:
      return Status("Model is proved either Infeasible or Unbounded", CPXMIP_INForUNBD);
    default:
      return Status("Unknown status code from Cplex");
    }
  }

  Status solve(IloNum * elapsed_time = NULL)
  {

    bool unsolved = false;
    try{
      if(!model_extracted)
	extractModel();
    
      if(elapsed_time != NULL)
	*elapsed_time = -solver.getImpl()->getCplexTime();

      if(! solver.solve() )
	{
	  *elapsed_time = 0;
	  unsolved = true;
	}

      if(elapsed_time != NULL)
	*elapsed_time += solver.getImpl()->getCplexTime();

    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    
    if (unsolved != true)
      model_solved = true;
    
    try
      {
	IloCplex::CplexStatus cplexStatus = solver.getCplexStatus();
	return handleCplexStatus(cplexStatus);
      }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    

  }

  Status readBasis(const char* filename)
  {
    try{
      solver.readBasis(filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status writeBasis(const char* filename)
  {
    try{
      solver.writeBasis(filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status readModel(const char* filename)
  {
    try{
      solver.importModel(model, filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status writeModel(const char* filename)
  {
    try{
      solver.exportModel(filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status readParam(const char* filename)
  {
    try{
      solver.readParam(filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status writeParam(const char* filename)
  {
    try{
      solver.writeParam(filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status readMipStart(const char* filename)
  {
    try{
      solver.readMIPStart(filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status writeMipStart(const char* filename)
  {
    try{
      solver.writeMIPStarts(filename,0,1);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  Status writeConflict(const char* filename)
  {
    try{
      solver.writeConflict(filename);
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
    return Status();
  }

  double getObjectiveValue()
  {
    if(!model_solved)
      return 0;
	    
    return solver.getObjValue();
  }

  Status getValues(NumericalArray& dest, const ExpressionArray& expr)
  {
    assert_equal(dest.shape(0), expr.shape(0));
    assert_equal(dest.shape(1), expr.shape(1));
	    
    // if(!model_solved)
    // 	return Status("Cannot get value; model not in a solved state.");
		    
    try{
      for(long i = 0; i < dest.shape(0); ++i)
	for(long j = 0; j < dest.shape(1); ++j)
	  dest(i,j) = solver.getValue(expr(i,j));
    }
    catch(IloException& e){
      return Status(e.getMessage());
    }
	    
    return Status();
  }

  long getNRows() const
  {
    if(solved())
      return solver.getNrows();
    else
      return 0;
  }

  long getNCols() const
  {
    if(solved())
      return solver.getNcols();
    else
      return 0;
  }

  long getNQCs() const
  {
    if(solved())
      return solver.getNQCs();
    else
      return 0;
  }

  long getNIterations() const
  {
    if(solved())
      return solver.getNiterations();
    else
      return 0;
  }

  double getBestObjValue() const
  {
    if(solved())
      return solver.getBestObjValue();
    else
      return 0;
  }

  double getCutoff() const
  {
    if(solved())
      return solver.getCutoff();
    else
      return 0;
  }

  double getMIPRelativeGap() const
  {
    if(solved())
      return solver.getMIPRelativeGap();
    else
      return 0;
  }

  int getNnodes() const
  {
    if(solved())
      return solver.getNnodes();
    else
      return 0;
  }


  bool solved() const { return model_solved; }
	
  string asString() const
  {
    ostringstream constraints;
    ostringstream objective;

    for(IloModel::Iterator it(model); it.ok(); ++it){
      if( (*it).isObjective())
	{
	  IloObjective obj = (*it).asObjective();
		    
	  objective << ( (obj.getSense() == IloObjective::Maximize) 
			 ? "maximize " : "minimize ")
		    << obj.getExpr() << " such that\n  ";
	}
      else
	constraints << (*it) << "\n  ";
    }
	    
    objective << constraints.str();
    return objective.str();
  }

private:

  void extractModel() 
  {
    env.setNormalizer(IloFalse);
    solver.extract(model);
    model_extracted = true;
    env.setNormalizer(IloTrue);
  }

  IloEnv env;
  IloModel model;
  IloCplex solver;
  IloObjective* current_objective;
  bool model_extracted;
  bool model_solved;
};

inline CPlexModelInterface::Status newCPlexModelInterface(CPlexModelInterface **cpx, IloEnv env)
{

  try {
    *cpx = new CPlexModelInterface(env);
    return CPlexModelInterface::Status();
  } catch(IloException& e) {
    return CPlexModelInterface::Status(e.getMessage());
  }
}

#endif /* _ARRAY_FUNCTIONS_H_ */



