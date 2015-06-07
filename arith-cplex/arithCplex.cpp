//
// arithCplex: A mixed integer programming model for generating mathematically valid
//             arithmetic equations over given constants and unknowns, respecting the
//             given hard and soft constraints; uses ILOG Concert library.
//
// Author: Ashish Sabharwal, May 2015
//
// Usage: arithCplex -h   ||   arithCplex --help
//
// The input is specified in a file formatted as follows:
//
//   constants             : 4 2
//   unknowns              : x
//   operators             : + - * / =
//   objtypes              : book person
//   constantOrUnknownType : 0 1 0
//   n                     : 5
//   answer                : 8
//
// Alternatively, a pre-generated model in .lp, .mps, or .sav file may be provided as input.
//

#include <algorithm>
#include <fstream>
#include <stack>
#include <vector>

#include <ilcplex/ilocplex.h>

ILOSTLBEGIN


// global constants;
const double epsilon = 0.000001;

// struct to capture various arithmetic expression formats
struct FormattedExpr {
  string postfix;
  string typedPostfix;
  string infix;
  string infixForSymPy;
  string answer;

  FormattedExpr() : answer("-") {}  // default numerical value of the unknown, i.e., the "answer"
};

// simple helper methods
inline bool isInt(const double x) {
  return fabs(x - round(x)) <= epsilon;
}

// forward declarations of methods
void parseCommandLine(int argc, char **argv);
void printCommandLine(int argc, char **argv);
void parseInputFile(const string infilename);
void setCplexParameters(IloCplex cplex);
void buildArithmeticModel(IloModel model, IloObjective obj, IloNumVarArray vars, IloRangeArray rngs);
void prettyPrintSoln(const IloIntArray intvals);
FormattedExpr getFormattedExpr(const IloIntArray intvals);
void solveExpressionsWithSymPy(vector<FormattedExpr> & expressions);


// model parameters to be read from the input file
vector<string> constants, unknowns, operators, objtypes;
vector<int>    constantOrUnknownType;
int            n;
double         trueAnswer;

// derived model parameters
int            l, k, p, q, m;
int            objtypeIdxNone = -1;   // index of "NONE" object type, if present
int            constantIdxZero = -1;  // index of "0" constant, if present
int            constantIdxOne = -1;   // index of "1" constant, if present
vector<double> constantValues;        // numerical values of the constants
vector<bool>   isIntConstant;         // indicator: is the i-th constant an integer?
bool           allIntConstants;       // indicator: are all constants integers?

// hardcoded model parameters
const int maxStackDepth = 10;

// (optional) model filename to import
char mip_filename[1024] = "";

// arithmetic question parameters
char arith_filename[1024] = "";  // file to import arithmetic model parameters from

// parameters defining what weight (possibly "infinite") to use with each arithmetic constraint
const int lowWt  = 1;
const int midWt  = 5;
const int highWt = 10;
const int infWt  = 1000000;

int wtNonOperatorsAtMostOnce           = infWt;
int wtStackDepthUpperBound             = infWt;
int wtExactlyOneUnknown                = infWt;
int wtNoTwoConsecutiveMultiplications  = infWt;
int wtNoTwoConsecutiveDivisions        = infWt;
int wtNoNegatives                      = infWt;

int wtTypeConsistency                  = highWt;
int wtEqualityFirstOrLast              = highWt;
int wtIntConstantsImplyIntUnknown      = highWt;

int wtPreserveOrderingInText           = midWt;

int wtUnknownFirstOrLast               = lowWt;
int wtEqualityNextToUnknown            = lowWt;

int wtHasAddition                      = 0; // lowWt;
int wtHasSubtraction                   = 0; // lowWt;
int wtHasMultiplication                = 0; // lowWt;
int wtHasDivision                      = 0; // lowWt;

// cplex parameters: defaults
int    param_nsolutions       = 25;
bool   param_printexpr        = true;
bool   param_printanswer      = true;
bool   param_printsoln        = false;
string param_savemodel        = "";
int    param_timelimit        = -1;
int    param_memlimit         = 2048;   // 2 GB by default
double param_epgap_percent    = 0;
int    param_threads          = 1;
int    param_mipinterval      = 0;  // let cplex decide logging interval

char   param_cplexrootalg[32] = "o";


// a cplex callback to print the current best found solution
ILOMIPINFOCALLBACK3(MIPInfoCallback, bool, isMinimization, IloCplex, cplex, IloNum, startTime) {
  if (!hasIncumbent())
    return;
  static double prevBest = (isMinimization ? 1e+100 : -1e+100);
  const IloNum currentBest = getIncumbentObjValue();
  if ((isMinimization && currentBest < prevBest) || (!isMinimization && currentBest > prevBest)) {
    prevBest = currentBest;
    cout << "New incumbent found at " << (cplex.getCplexTime() - startTime) << " sec : " << currentBest << endl;
  }
}


//////////////////////////////////////////////////
int main (int argc, char **argv) {

  IloEnv env;

  try {
    IloCplex cplex(env);
    IloTimer timer(env);

    const IloNum startTime = cplex.getCplexTime();

    parseCommandLine(argc, argv);
    printCommandLine(argc, argv);

    cout << "Starting IloTimer" << endl;
    timer.start();

    // parse arithmetic model parameters
    if (strlen(arith_filename) > 0)
      parseInputFile(arith_filename);

    // set cplex paramters
    setCplexParameters(cplex);

    // import or build the model
    IloModel       model(env);
    IloObjective   obj(env);
    IloNumVarArray corevars(env);
    IloRangeArray  rngs(env);

    if (strlen(mip_filename) > 0) {
      cout << "------- Importing the model -------" << endl;
      cplex.importModel(model, mip_filename, obj, corevars, rngs);
      cout << "Model imported at " << (cplex.getCplexTime() - startTime) << " sec" << endl;
    }
    else {
      cout << "------- Building the model -------" << endl;
      buildArithmeticModel(model, obj, corevars, rngs);
    }

    int nvars = corevars.getSize();
    cout << "Number of Core Variables: " << nvars << endl;

    cplex.extract(model);
    cout << "Model extracted at " << (cplex.getCplexTime() - startTime) << " sec" << endl;

    // save the MIP model to a file, if desired
    if (!param_savemodel.empty()) {
      cout << endl << "Saving generated MIP model to " << param_savemodel << endl << endl;
      cplex.exportModel(param_savemodel.c_str());
    }

    // find out whether it is a minimization problem or a maximization one
    bool isMinimization = true;
    if (cplex.getObjective().getSense() == IloObjective::Maximize)
      isMinimization = false;

    // ask cplex to use MIPInfoCallback
    cplex.use(MIPInfoCallback(env, isMinimization, cplex, startTime));

    cout << "------- Solving the extracted model -------" << endl;
    //const bool solutionFound = cplex.solve();
    const bool solutionFound = (param_nsolutions > 1 ? cplex.populate() : cplex.solve());
    const int nSolutionsFound = cplex.getSolnPoolNsolns();

    cout << "Stopping IloTimer" << endl;
    timer.stop();
    const IloNum endTime = cplex.getCplexTime();

    cout << "-------------------------------------------" << endl;
    printCommandLine(argc, argv);
    cout << "-------------------------------------------" << endl;
    cout << "Number of cplex nodes        = " << cplex.getNnodes() << endl;
    cout << "Number of cplex iterations   = " << cplex.getNiterations() << endl;
    cout << "Aggregated CPU time          = " << timer.getTime() << " seconds" << endl;
    cout << "Elapsed wall clock time      = " << endTime - startTime << " seconds" << endl;
    cout << "Number of threads used       = " << param_threads << endl;
    cout << "Solution status              = " << cplex.getStatus() << endl;

    if (solutionFound) {
      cout << "Solution value               = " << cplex.getObjValue() << endl;
      cout << "Optimality Gap (in %)        = " << fabs((cplex.getBestObjValue() - cplex.getObjValue()) / (1.0 * cplex.getObjValue())) * 100 << endl;
      //cout << "Maximum bound violation      = " << cplex.getQuality(IloCplex::MaxPrimalInfeas) << endl;
    }
    cout << endl << "parameters: n=" << n << " l=" << l << " k=" << k << " p=" << p << " q=" << q << " m=" << m << endl;
    if (solutionFound) {
      cout << "TOTAL " << nSolutionsFound << " solutions found" << endl;
      int nAllowedSolutionsFound = 0;
      if (param_printexpr || param_printanswer || param_printsoln) {
        IloNumArray objValues(env, nSolutionsFound);
        vector<pair<IloNum,unsigned> > sortedIndex(nSolutionsFound);  // to sort solutions by objective value
        vector<FormattedExpr> expressions(nSolutionsFound);
        // extract all solutions as formatted expressions
        for (int s=0; s<nSolutionsFound; s++) {
          IloNumArray vals(env);
          cplex.getValues(vals, corevars, s);
          // convert solution values to integers; note: simple int cast may lead to errors!
          IloIntArray intvals(env, vals.getSize());
          for (int i=0; i<vals.getSize(); i++)
            intvals[i] = IloRound(vals[i]);  // use IloRound rather than std::round
          if (param_printsoln)
            prettyPrintSoln(intvals);
          objValues[s] = cplex.getObjValue(s);
          expressions[s] = getFormattedExpr(intvals);
          sortedIndex[s] = pair<IloNum,unsigned> (objValues[s],s);
        }
        // sort solutions by increasing objective value
        std::stable_sort(sortedIndex.begin(), sortedIndex.end());
        // evaluate all expressions with a single call to Python's SymPy package
        if (param_printanswer)
          solveExpressionsWithSymPy(expressions);
        // print expressions if desired, in sorted order
        if (param_printexpr) {
          cout << "SOLN: CORRECT | POS/NEG | INT/FRA | OBJ-SCORE | TRUE-ANS | ANS | INFIX | POSTFIX | TYPED-POSTFIX" << endl;
          for (int s=0; s<nSolutionsFound; s++) {
            const int sId = sortedIndex[s].second;
            const FormattedExpr & expr = expressions[sId];
            const double answerValue = atof(expr.answer.c_str());
            const bool isCorrect = (fabs(answerValue - trueAnswer) < epsilon);
            const bool isAnswerNegative = answerValue < 0;
            const bool isAnswerInteger = isInt(answerValue);
            if (!isAnswerNegative && (!allIntConstants || wtIntConstantsImplyIntUnknown == 0 || isAnswerInteger)) {
              ++nAllowedSolutionsFound;
              cout << "EXPR: " << isCorrect
                   << " | " << (isAnswerNegative ? "NEG" : "POS")
                   << " | " << (isAnswerInteger ? "INT" : "FRA")
                   << " | " << objValues[sId]
                   << " | " << trueAnswer
                   << " | " << expr.answer
                   << " | " << expr.infix
                   << " | " << expr.postfix
                   << " | " << expr.typedPostfix
                   << endl;
            }
          }
        }
      }
      const string solnProperty = string(" non-negative")
        + (allIntConstants && wtIntConstantsImplyIntUnknown != 0 ? ", integer-valued " : " ");
      cout << "NET " << nAllowedSolutionsFound << solnProperty << "solutions found out of "
           << nSolutionsFound << " total solutions" << endl;
    }

    cout << "-------------------------------------------" << endl;
    cout << "RESULT:"
         << " NODES " << cplex.getNnodes()
         << " | ITERATIONS " << cplex.getNiterations()
         << " | CPUTIME " << timer.getTime()
         << " | WALLTIME " << endTime - startTime
         << " | THREADS " << param_threads
         << " | STATUS " << cplex.getStatus();
    if (solutionFound)
      cout << " | SOLUTION " << cplex.getObjValue()
           << " | OPTGAP " << fabs((cplex.getBestObjValue() - cplex.getObjValue()) / (1.0 * cplex.getObjValue())) * 100;
    else
      cout << " | SOLUTION - | OPTGAP -";
    cout << endl;

    //try {     // basis may not exist
    //  IloCplex::BasisStatusArray cstat(env);
    //  cplex.getBasisStatuses(cstat, vars);
    //  cout << "Basis statuses               = " << cstat << endl;
    //}
    //catch (...) {
    //}


  }
  catch (IloException& e) {
    cerr << "Concert exception caught: " << e << endl;
  }
  catch (...) {
    cerr << "Unknown exception caught" << endl;
  }

  env.end();
  return 0;
}  // END main


// print usage
void printUsage(const char *progname, ostream & str = cout) {
   str << endl
       << "USAGE: " << progname << " [options] arith-filename" << endl
       << endl
       << " where arith-filename is a file with parameters for building an arithmetic problem model" << endl
       << " and options are:" << endl << endl
       << "   -h, --help        print this usage" << endl
       << endl
       << "   --model file      an ILP model file with extension MPS, SAV, or LP (lower case ok; .gz ok; default: none)" << endl
       << endl
       << "   --timelimit num   time limit in seconds (elapsed/wall time; default: none)" << endl
       << "   -t num            same as --timelimit num" << endl
       << "   --memlimit num    memory limit in MB (default: 2048)" << endl
       << "   -m num            same as --memlimit num" << endl
       << "   --threads num     number of parallel threads to use (default: 1)" << endl
       << "   -g percent        stop optimization when optimality gap reaches provided value" << endl
       << "   --log num         report log after every num nodes (default: 0, cplex decides)" << endl
       << endl
       << "   -s num            same as --solutions num" << endl
       << "   --solutions num   number of solutions to find (default: 25)" << endl
       << "   --noprintexpr     do not print arithmetic expressions found (forces --noprintexpr; default: on)" << endl
       << "   --noprintanswer   do not print answer to arithmetic problem (default: on)" << endl
       << "   --printsoln       print solution (default: off)" << endl
       << "   --savemodel file  save generated MIP model to file (default: none)" << endl
       << endl
       << "   --rootalg alg     algorithm to use for solving LP" << endl
       << "                        o default, p primal simplex, d dual simplex" << endl
       << "                        b barrier, h barrier with crossover," << endl
       << "                        n network simplex, s sifting, c concurrent" << endl
       << "   -a alg            same as --rootalg alg" << endl
       << endl;
} // END printUsage

// parse command line arguments
void parseCommandLine(int argc, char **argv) {
  bool flag_read_filename = false;
  for (int i=1; i<argc; i++) {
    char *optionName = argv[i];
    if (!strcmp(optionName, "--help") || !strcmp(optionName, "-h")) {
      printUsage(argv[0], cout);
      exit(0);
    }
    else if (!strcmp(optionName, "--model"))
      strcpy(mip_filename, argv[++i]);
    else if (!strcmp(optionName, "--log"))
      param_mipinterval = atoi(argv[++i]);
    else if (!strcmp(optionName, "--rootalg") || !strcmp(optionName, "-a"))
      strcpy(param_cplexrootalg, argv[++i]);
    else if (!strcmp(optionName, "--solutions") || !strcmp(optionName, "-s"))
      param_nsolutions = atoi(argv[++i]);
    else if (!strcmp(optionName, "--noprintexpr")) {
      param_printexpr = false;
      param_printanswer = false;   // force param_printanswer to be false
    }
    else if (!strcmp(optionName, "--noprintanswer"))
      param_printanswer = false;
    else if (!strcmp(optionName, "--printsoln"))
      param_printsoln = true;
    else if (!strcmp(optionName, "--savemodel"))
      param_savemodel = argv[++i];
    else if (!strcmp(optionName, "--timelimit") || !strcmp(optionName, "-t"))
      param_timelimit = atoi(argv[++i]);
    else if (!strcmp(optionName, "--memlimit") || !strcmp(optionName, "-m"))
      param_memlimit = atoi(argv[++i]);
    else if (!strcmp(optionName, "--threads"))
      param_threads = atoi(argv[++i]);
    else if (!strcmp(optionName, "-g"))
      param_epgap_percent = atof(argv[++i]);
    else if (optionName[0] == '-') {
      cerr << endl << "ERROR: invalid command-line option \"" << optionName << '\"' << endl;
      printUsage(argv[0], cerr);
      exit(-1);
    }
    else {
      // this must be the input filename
      if (flag_read_filename == true) {
        printUsage(argv[0], cerr);
        exit(-1);
      }
      strcpy(arith_filename, optionName);
      flag_read_filename = true;
    }
  }
  if (flag_read_filename == false) {
    cout << "No arithmetic input file specified" << endl;
    printUsage(argv[0], cerr);
    exit(1);
  }
}

void printCommandLine(int argc, char **argv) {
  cout << "Command line: ";
  for (int i=0; i<argc-1; i++)
    cout << argv[i] << ' ';
  cout << argv[argc-1] << endl;
}

void parseInputFile(const string infilename) {
  cout << "Reading input from file " << infilename << endl;
  std::ifstream infile(infilename);
  string line;
  while (std::getline(infile, line)) {
    cout << line << endl;
    if (line.empty() || line.substr(0, 2) == "//")
      continue;
    std::istringstream iss(line);
    string header, word;
    iss >> header;
    if (header == "constants") {
      iss >> word; assert(word == ":");
      while (iss >> word) {
        // drop commas from the word, e.g., 2,715 turns into 2715
        size_t pos = 0;
        while ((pos = word.find(',', pos)) != string::npos)
          word.erase(pos, 1);
        if (word == "0")
          constantIdxZero = constants.size();
        else if (word == "1")
          constantIdxOne = constants.size();
        constants.push_back(word);
        const double numericalValue = atof(word.c_str());
        constantValues.push_back(numericalValue);
        isIntConstant.push_back(isInt(numericalValue));
      }
    }
    else if (header == "unknowns") {
      iss >> word; assert(word == ":");
      while (iss >> word) unknowns.push_back(word);
    }
    else if (header == "operators") {
      iss >> word; assert(word == ":");
      while (iss >> word) operators.push_back(word);
      assert(operators.size() == 5);
    }
    else if (header == "objtypes") {
      iss >> word; assert(word == ":");
      while (iss >> word) {
        if (word == "NONE")
          objtypeIdxNone = objtypes.size();
        objtypes.push_back(word);
      }
    }
    else if (header == "constantOrUnknownType") {
      iss >> word; assert(word == ":");
      while (iss >> word) constantOrUnknownType.push_back(atoi(word.c_str()));
    }
    else if (header == "n") {
      iss >> word; assert(word == ":");
      iss >> word;
      n = atoi(word.c_str());
    }
    else if (header == "answer") {
      iss >> word; assert(word == ":");
      iss >> word;
      trueAnswer = atof(word.c_str());
    }
    else {
      cerr << "ERROR: unrecognized header in line: " << line << endl;
      exit(1);
    }
  }
  infile.close();
  cout << endl;

  // set values of dependent model parameters
  l = constants.size();
  k = unknowns.size();
  p = operators.size();
  q = objtypes.size();
  m = l + k + p;

  // check whether all specified constants are integers
  allIntConstants = (std::find(isIntConstant.begin(), isIntConstant.end(), false) == isIntConstant.end());
  if (allIntConstants)
    cout << "All parsed constants are integers. Forcing answer to be an integer." << endl << endl;

  // special treatment for objtype "dozen"
  for (unsigned i=0; i<objtypes.size(); i++) {
    string str = objtypes[i].substr(0,5);
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    if (str == "dozen") {
      cout << "Detected object type starting with 'dozen'." << endl
           << "  Turning off NoTwoConsecutiveMultiplications constraint." << endl << endl;
      wtNoTwoConsecutiveMultiplications = 0;
      break;
    }
  }
}

void setCplexParameters(IloCplex cplex) {
  // main search algorithm
  if (!strcmp(param_cplexrootalg, "o")) cplex.setParam(IloCplex::RootAlg, IloCplex::AutoAlg);
  else if (!strcmp(param_cplexrootalg, "p")) cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
  else if (!strcmp(param_cplexrootalg, "d")) cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
  else if (!strcmp(param_cplexrootalg, "b")) {
    cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier);
    cplex.setParam(IloCplex::BarCrossAlg, IloCplex::NoAlg);
  }
  else if (!strcmp(param_cplexrootalg, "h")) cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier);
  else if (!strcmp(param_cplexrootalg, "n")) cplex.setParam(IloCplex::RootAlg, IloCplex::Network);
  else if (!strcmp(param_cplexrootalg, "s")) cplex.setParam(IloCplex::RootAlg, IloCplex::Sifting);
  else if (!strcmp(param_cplexrootalg, "c")) cplex.setParam(IloCplex::RootAlg, IloCplex::Concurrent);
  else {
    cerr << "ERROR: unrecognized value of option -a / --rootalg" << endl;
    throw(-1);
  }

  // search limits
  if (param_timelimit > 0) {
    cplex.setParam(IloCplex::TiLim, param_timelimit);
    cplex.setParam(IloCplex::ClockType, 2);  // 0: auto; 1: cpu; 2: elapsed/wall (default)
  }

  // parallel mode
  cplex.setParam(IloCplex::Threads, param_threads);
  if (param_threads > 1)
    cplex.setParam(IloCplex::ParallelMode, 1); // 1: deterministic, -1: opportunistic; 0: auto

  // output log frequency
  cplex.setParam(IloCplex::MIPInterval, param_mipinterval);

  // memory related parameters
  cplex.setParam(IloCplex::NodeFileInd, 2); // write to disk but don't compress
  cplex.setParam(IloCplex::WorkMem, param_memlimit);

  // solution quality
  if (param_epgap_percent > 0)
    cplex.setParam(IloCplex::EpGap, param_epgap_percent / 100.0);

  // multiple solutions
  if (param_nsolutions > 1) {
    cplex.setParam(IloCplex::PopulateLim, 10000);  // explore at most so many solutions
    cplex.setParam(IloCplex::SolnPoolCapacity, param_nsolutions);  // keep at most so many in the pool
    cplex.setParam(IloCplex::SolnPoolReplace, 1);  // discard solutions based on obj value
    cplex.setParam(IloCplex::SolnPoolIntensity, 4);  // enumerate all
  }
}

void buildArithmeticModel(IloModel model, IloObjective obj, IloNumVarArray vars, IloRangeArray rngs) {
  const IloEnv env = model.getEnv();

  // add placeholder variables and objective to the model
  IloIntVarArray intVars(env);
  model.add(intVars);

  IloExpr objective(env);  // to be populated, then set as the expression for obj
  obj.setSense(IloObjective::Minimize);
  model.add(obj);

  // create model variables
  IloIntVarArray x(env, n, 0, m-1);      x.setNames("x");
  IloIntVarArray c(env, n, 0, 1);        c.setNames("c");
  IloIntVarArray u(env, n, 0, 1);        u.setNames("u");
  IloIntVarArray o(env, n, 0, 1);        o.setNames("o");
  IloIntVarArray t(env, n, 0, q-1);      t.setNames("t");
  IloIntVarArray d(env, n, 0, n-1);      d.setNames("d");
  IloIntVarArray op1idx(env, n, 0, n-1); op1idx.setNames("op1idx");
  IloIntVarArray op1x(env, n, 0, m-1);   op1x.setNames("op1x");
  IloIntVarArray op1u(env, n, 0, 1);     op1u.setNames("op1u");
  IloIntVarArray op1o(env, n, 0, 1);     op1o.setNames("op1o");
  IloIntVarArray op1t(env, n, 0, q-1);   op1t.setNames("op1t");
  intVars.add(x);
  intVars.add(c);
  intVars.add(u);
  intVars.add(o);
  intVars.add(t);
  intVars.add(d);
  intVars.add(op1idx);
  intVars.add(op1x);
  intVars.add(op1u);
  intVars.add(op1o);
  intVars.add(op1t);

  // store core variables to vars array to be returned
  vars.add(intVars.toNumVarArray());

  // add constraints

  // definition of C variables
  for (int i=0; i<n; i++)
    //model.add(x[i] + l*c[i] >= l);
    //model.add(x[i] + (m-l)*c[i] <= m-1);
    model.add((c[i] == 1) == (x[i] <= l-1));

  // definition of U variables (indirect)
  for (int i=0; i<n; i++)
    model.add(c[i] + u[i] + o[i] == 1);

  // definition of O variables
  for (int i=0; i<n; i++)
    //model.add(p*o[i] >= x[i] - l - k + 1);
    //model.add((l+k)*o[i] <= x[i]);
    model.add((o[i] == 1) == (x[i] >= l+k));

  // definition of D variables
  model.add(d[0] == 0);
  for (int i=1; i<n; i++)
    model.add(d[i] == d[i-1] - 2*o[i] + 1);

  // definition of OP1 variables
  // op1idx[i] = argmax_{j<i-1} {d[j] == d[i]}
  model.add(op1idx[0] == 0);  // irrelevant
  model.add(op1idx[1] == 0);  // irrelevant
  for (int i=2; i<n; i++) {
    model.add(op1idx[i] <= i-2);
    model.add(op1idx[i] <= (n-1)*o[i]);  // irrlevant if o[i] == 0
    for (int j=0; j<i-1; j++) {
      model.add((d[i] == d[j]) <= (1-o[i]) + (op1idx[i] >= j));
      model.add((d[i] != d[j]) <= (1-o[i]) + (op1idx[i] != j));
    }
  }
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      model.add((op1idx[i] == j) <= (op1x[i] == x[j]));
      model.add((op1idx[i] == j) <= (op1u[i] == u[j]));
      model.add((op1idx[i] == j) <= (op1o[i] == o[j]));
      model.add((op1idx[i] == j) <= (op1t[i] == t[j]));
    }
  }

  // validity of the postfix expression: stack depth
  model.add(o[0] == 0);
  //for (int i=1; i<n-1; i++) model.add(d[i] >= 0);  // implied by the domain of D vars
  model.add(d[n-1] == 0);
  // validity: equality must occur exactly once and must be the "top" level operator,
  //  i.e., at the very end in postfix
  model.add(x[n-1] == m-1);
  for (int i=0; i<n-1; i++)
    model.add(x[i] <= m-2);

  // simplicity: "0" is not an operand of + or -, or the second operand of /
  for (int i=2; i<n-1; i++) {
    model.add((x[i] == l+k) <= (op1x[i] != constantIdxZero));
    model.add((x[i] == l+k) <= (x[i-1] != constantIdxZero));
    model.add((x[i] == l+k+1) <= (op1x[i] != constantIdxZero));
    model.add((x[i] == l+k+1) <= (x[i-1] != constantIdxZero));
    model.add((x[i] == l+k+3) <= (x[i-1] != constantIdxZero));
  }

  // simplicity: "1" is not an operand of * or the second operand of /
  for (int i=2; i<n-1; i++) {
    //model.add((x[i] == l+k+2) <= (op1x[i] != constantIdxOne));
    //model.add((x[i] == l+k+2) <= (x[i-1] != constantIdxOne));
    model.add((x[i] == l+k+3) <= (x[i-1] != constantIdxOne));
  }

  // symmetry breaking for commutative operators:
  // a + preceded by two constants must have those constants in the given order
  for (int i=2; i<n-1; i++)
    model.add((x[i] == l+k) + c[i-1] + c[i-2] <= 2 + (x[i-2] <= x[i-1] - 1));
  // a * preceded by two constants must have those constants in the given order
  for (int i=2; i<n-1; i++)
    model.add((x[i] == l+k+2) + c[i-1] + c[i-2] <= 2 + (x[i-2] <= x[i-1] - 1));
  // disable expressions of the form LHS = x; these are equivalent to x = RHS
  model.add(u[n-2] == 0);

  // (optional) TypeConsistency
  if (wtTypeConsistency == 0) {
    for (int i=0; i<n; i++)
      model.add(t[i] == 0);
  }
  else {
    // enforce (potentially error prone) types of base entities
    for (int i=0; i<n; i++) {
      for (int j=0; j<l+k; j++) {
        if (constantOrUnknownType[j] != objtypeIdxNone) {  // type NONE is a wildcard
          IloNumVar slack(env);
          model.add((x[i] == j) <= (t[i] == constantOrUnknownType[j]) + slack);
          objective += slack * wtTypeConsistency;
        }
      }
    }
    // define semantics of operators;
    // if x[i] is an operator, its operands are op1idx[i] and x[i-1]
    for (int i=2; i<n; i++) {
      const IloIntVar z(env, -l-k, m-1-l-k);
      model.add(z == x[i]-l-k);
      model.add((z == 0) <= (t[i-1] == op1t[i]));
      model.add((z == 1) <= (t[i-1] == op1t[i]));
      model.add((z == 4) <= (t[i-1] == op1t[i]));
      model.add((z == 0) <= (t[i] == t[i-1]));
      model.add((z == 1) <= (t[i] == t[i-1]));
      model.add((z == 4) <= (t[i] == t[i-1]));
      model.add((z == 2) <= (t[i-1] != op1t[i]));
      model.add((z == 2) <= (t[i] == t[i-1]) + (t[i] == op1t[i]));
      model.add((z == 3) <= (t[i-1] != op1t[i]));
      model.add((z == 3) <= (t[i] == op1t[i]));
    }
  }

  // (optional) NonOperatorsAtMostOnce
  if (wtNonOperatorsAtMostOnce != 0) {
    if (wtPreserveOrderingInText != infWt) {   // otherwise unnecessary
      for (int i=0; i<n-2; i++) {
        for (int j=i+1; j<n-1; j++) {
          IloNumVar slack(env);
          model.add((x[i] == x[j]) <= o[i] + slack);
          objective += slack * wtNonOperatorsAtMostOnce;
        }
      }
    }
  }

  // (optional) simplicity: StackDepthUpperBound
  if (wtStackDepthUpperBound != 0) {
    for (int i=1; i<n-1; i++) {
      IloNumVar slack(env);
      model.add(d[i] <= maxStackDepth - 1 + slack);
      objective += slack * wtStackDepthUpperBound;
    }
  }

  // (optional) simplicity: EqualityFirstOrLast
  // In postfix, this means at least one of the operands of equality (which appears
  // always at position n-1 due to validity requirements) must not be an operator
  if (wtEqualityFirstOrLast != 0) {
    IloNumVar slack(env);
    model.add(op1o[n-1] + o[n-2] <= 1 + slack);
    objective += slack * wtEqualityFirstOrLast;
  }

  // (optional) simplicity: NoTwoConsecutiveMultiplications
  if (wtNoTwoConsecutiveMultiplications != 0) {
    for (int i=2; i<n-2; i++) {
      IloNumVar slack(env);
      model.add((x[i] == l+k+2) + (x[i+1] == l+k+2) + (x[i+2] == l+k+2) <= 1 + slack);
      objective += slack * wtNoTwoConsecutiveMultiplications;
    }
  }

  // (optional) simplicity: NoTwoConsecutiveDivisions
  if (wtNoTwoConsecutiveDivisions != 0) {
    for (int i=2; i<n-2; i++) {
      IloNumVar slack(env);
      model.add((x[i] == l+k+3) + (x[i+1] == l+k+3) + (x[i+2] == l+k+3) <= 1 + slack);
      objective += slack * wtNoTwoConsecutiveDivisions;
    }
  }

  // (optional) simplicity: NoNegatives
  if (wtNoNegatives != 0) {
    // TO-DO!!!
  }

  // (optional) simplicity: IntConstantsImplyIntUnknowns
  if (wtIntConstantsImplyIntUnknown != 0) {
    // TO-DO!!!
  }

  // (optional) PreserveOrderingInText
  if (wtPreserveOrderingInText != 0) {
    // apply penalty at the entity level, not entity-pair level (otherwise it may
    // grow too quickly)
    for (int i=0; i<n-2; i++) {
      IloNumVar slack(env);
      for (int j=i+1; j<n-1; j++) {
        model.add(x[i] <= (x[j]-1) + (m-1)*(1-c[i]) + slack);
      }
      objective += slack * wtPreserveOrderingInText;
    }
  }

  // (optional) ExactlyOneUnknown; note: at least one unknown is a must; at most once is optional
  model.add(IloSum(u) >= 1);
  if (wtExactlyOneUnknown != 0) {
    IloNumVar slack(env);
    model.add(IloSum(u) <= 1 + slack);
    objective += slack * wtExactlyOneUnknown;
  }

  // (optional) UnknownFirstOrLast
  if (wtUnknownFirstOrLast != 0) {
    IloIntVar y0(env, 0, 1);
    IloIntVar y1(env, 0, 1);
    model.add(y0 <= u[0]);
    for (int i=1; i<n-1; i++)
      model.add(y0 <= d[i]);
    model.add(y1 <= u[n-2]);
    IloNumVar slack(env);
    model.add(y0 + y1 + slack >= 1);
    objective += slack * wtUnknownFirstOrLast;
  }

  // (optional) EqualityNextToUnknown
  if (wtEqualityNextToUnknown != 0) {
    for (int i=2; i<n; i++) {
      IloNumVar slack(env);
      model.add(x[i] - m + 2 <= u[i-1] + op1u[i] + slack);
      objective += slack * wtEqualityNextToUnknown;
    }
  }

  // (optional) HasAddition, HasSubtraction, HasMultiplication, HasDivision
  if (wtHasAddition != 0) {
    IloExpr y(env);
    for (int i=2; i<n; i++)
      y += (x[i] == l+k);
    IloNumVar slack(env);
    model.add(y + slack >= 1);
    objective += slack * wtHasAddition;
  }
  if (wtHasSubtraction != 0) {
    IloExpr y(env);
    for (int i=2; i<n; i++)
      y += (x[i] == l+k+1);
    IloNumVar slack(env);
    model.add(y + slack >= 1);
    objective += slack * wtHasSubtraction;
  }
  if (wtHasMultiplication != 0) {
    IloExpr y(env);
    for (int i=2; i<n; i++)
      y += (x[i] == l+k+2);
    IloNumVar slack(env);
    model.add(y + slack >= 1);
    objective += slack * wtHasMultiplication;
  }
  if (wtHasDivision != 0) {
    IloExpr y(env);
    for (int i=2; i<n; i++)
      y += (x[i] == l+k+3);
    IloNumVar slack(env);
    model.add(y + slack >= 1);
    objective += slack * wtHasDivision;
  }

  // set the expression to be optimized
  model.add(objective <= infWt - 1);  // should not violate hard constraints
  obj.setExpr(objective);
}

void prettyPrintSoln(const IloIntArray intvals) {
  int i=0;
  cout << endl << "x: ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "c: ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "u: ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "o: ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "t: ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "d: ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "op1idx: ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "op1x  : ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "op1u  : ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "op1o  : ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl << "op1t  : ";
  for (int l=0; l<n; l++) cout << intvals[i++] << ' ';
  cout << endl;
}

FormattedExpr getFormattedExpr(const IloIntArray intvals) {
  FormattedExpr expr;
  std::stack<string> stk;
  for (int i=0; i<n; i++) {
    const int intval = intvals[i];
    string symbol;
    if (intval < l) {
      symbol = constants[intval];
      stk.push(symbol);
    }
    else if (intval < l+k) {
      symbol = unknowns[intval-l];
      stk.push(symbol);
    }
    else if (intval < l+k+m) {
      symbol = operators[intval-l-k];
      assert(stk.size() >= 2);
      const string v2 = stk.top(); stk.pop();
      const string v1 = stk.top(); stk.pop();
      if (i < n-1)
        stk.push('(' + v1 + symbol + v2 + ')');
      else {
        assert(symbol == "=");
        assert(stk.empty());
        expr.infix = v1 + symbol + v2;
        expr.infixForSymPy = v1 + '-' + v2;
      }
    }
    else {
      cerr << "ERROR: unexpected variable value " << intval << " (out of domain)." << endl;
      exit(1);
    }
    expr.postfix += symbol + ' ';
    const int objtype = intvals[4*n + i];
    expr.typedPostfix += symbol + ':' + objtypes[objtype] + ' ';
  }
  expr.postfix.erase(expr.postfix.size()-1);
  expr.typedPostfix.erase(expr.typedPostfix.size()-1);
  return expr;
}

void solveExpressionsWithSymPy(vector<FormattedExpr> & expressions) {
  if (expressions.empty())
    return;
  // create a comma-separated list of single-quoted expressions to feed to SymPy
  // e.g. '(4+2)-x','x-2'
  string exprList = "'" + expressions[0].infixForSymPy + "'";
  for (unsigned i=1; i<expressions.size(); i++)
    exprList += ",'" + expressions[i].infixForSymPy + "'";
  // invoke SymPy through a pipe and capture the output
  const string pythonCmd = string("python -c \"")
    + "import sympy\n"
    + "x=sympy.Symbol('x')\n"
    + "def mysolve(e):\n"
    + "  try:\n"
    + "    return str((sympy.solve(e+'*1.0'))[0])\n"
    + "  except:\n"
    + "    return '-'\n"
    + "print [mysolve(e).rstrip('0').rstrip('.') for e in [" + exprList + "]]\"";
  //cout << pythonCmd << endl;
  FILE* pipe = popen(pythonCmd.c_str(), "r");
  string result;
  char buffer[128];
  while (!feof(pipe))
    if (fgets(buffer, 128, pipe) != NULL)
      result += buffer;
  pclose(pipe);
  // parse output into answers
  // result should be in the format [a1, a2, ..., ak]
  assert(result.size() >= 3);
  result.erase(0, 1);   // drop leading '['
  result.erase(result.size() - 2);   // drop trailing ']' and newline
  stringstream ss(result);
  string answer;
  int i=0;
  while (ss.good()) {
    getline(ss, answer, ',');
    if (answer[0] == ' ')
      answer.erase(0, 1);
    answer.erase(0, 1);  // drop leading '
    answer.erase(answer.size()-1);  // drop trailing '
    expressions[i].answer = answer;
    ++i;
  }
  assert(i == expressions.size());
}

