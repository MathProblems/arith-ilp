# arith-cplex
CPLEX based Integer Linear Programming models for arithmetic problems

Compilation
-----------

This project requires `CPLEX` libraries and optionally also uses the `SymPy` package in `Python` for evaluating expressions (can be turned off with `--noprintanswer`). To compile:
 
* copy `Makefile-SAMPLE` to your own local `Makefile` (do not check this
  in, as this will be computer specific)

* edit your local `Makefile` and set `CPLEXDIR` and `CONCERTDIR` variables appropriately

* run `make`

This should produce an executable named `arithCplex`.


Usage
-----

**Running the CPLEX model on one input**

Type `arithCplex -h` for a list of all options and their default values. Options include multi-threaded runs (e.g., `--threads 8`), the desired number of solutions (e.g., `-s 100`), disabling expression evaluation using `Python`, etc.

Basic usage: `arithCplex inputfile` where `inputfile` is a plain text file containing parameters to build an arithmetic expression model, such as `example.txt` or the following:

```
constants : 8 3
unknowns : x
operators : + - * / =
objtypes : calorie bar
constantOrUnknownType : 0 1 0
n : 5
answer : 24
```

Note that the syntax of this input file is rather brittle. It uses white-space as delimiters and, in particular, expects a white-space on both sides of the `:` character.

The output includes a number of CPLEX related info as well as the number of solutions found. For each solution, it prints information in the format. For the included `example.txt` arithmetic problem file, the output looks like:

```
parameters: n=5 l=2 k=1 p=5 q=1 m=8
TOTAL 22 solutions found
SOLN: CORRECT | POS/NEG | INT/FRA | OBJ-SCORE | TRUE-ANS | ANS | INFIX | POSTFIX | TYPED-POSTFIX
EXPR: 0 | POS | INT | 0 | 43 | 97 | x=(70+27) | x 70 27 + = | x:seashell 70:seashell 27:seashell +:seashell =:seashell
EXPR: 1 | POS | INT | 0 | 43 | 43 | x=(70-27) | x 70 27 - = | x:seashell 70:seashell 27:seashell -:seashell =:seashell
EXPR: 1 | POS | INT | 1 | 43 | 43 | x=(70-27) | x 70 27 - = | x:seashell 70:seashell 27:seashell -:seashell =:seashell
EXPR: 0 | POS | INT | 1 | 43 | 97 | x=(70+27) | x 70 27 + = | x:seashell 70:seashell 27:seashell +:seashell =:seashell
EXPR: 0 | POS | INT | 2 | 43 | 97 | 70=(x-27) | 70 x 27 - = | 70:seashell x:seashell 27:seashell -:seashell =:seashell
EXPR: 0 | POS | INT | 2 | 43 | 97 | (x-70)=27 | x 70 - 27 = | x:seashell 70:seashell -:seashell 27:seashell =:seashell
EXPR: 1 | POS | INT | 2 | 43 | 43 | 70=(27+x) | 70 27 x + = | 70:seashell 27:seashell x:seashell +:seashell =:seashell
EXPR: 1 | POS | INT | 2 | 43 | 43 | 70=(x+27) | 70 x 27 + = | 70:seashell x:seashell 27:seashell +:seashell =:seashell
EXPR: 1 | POS | INT | 2 | 43 | 43 | (70-x)=27 | 70 x - 27 = | 70:seashell x:seashell -:seashell 27:seashell =:seashell
EXPR: 0 | POS | INT | 12 | 43 | 97 | 27=(x-70) | 27 x 70 - = | 27:seashell x:seashell 70:seashell -:seashell =:seashell
EXPR: 1 | POS | INT | 12 | 43 | 43 | (27+x)=70 | 27 x + 70 = | 27:seashell x:seashell +:seashell 70:seashell =:seashell
EXPR: 0 | POS | INT | 12 | 43 | 97 | (x-27)=70 | x 27 - 70 = | x:seashell 27:seashell -:seashell 70:seashell =:seashell
EXPR: 1 | POS | INT | 12 | 43 | 43 | (x+27)=70 | x 27 + 70 = | x:seashell 27:seashell +:seashell 70:seashell =:seashell
EXPR: 1 | POS | INT | 12 | 43 | 43 | 27=(70-x) | 27 70 x - = | 27:seashell 70:seashell x:seashell -:seashell =:seashell
NET 14 non-negative, integer-valued solutions found out of 22 total solutions
```


**Running the CPLEX model on all `*.txt` files in a folder**

For experimentation, the included `runall.sh` script in the `scripts` folder may be used to run the CPLEX model on all `*.txt` files included in a folder. The output for each arithmetic problem file `file.txt` is then captured in `file.txt.out`.

Usage: `./runall.sh datafolder`


**Splitting an joint input file for a dataset into multiple `*.txt` files for `arithCplex`**

The `splitQuestions.sh` script in the `scripts` folder may be used to split a joint file, such as `emnlp.curated.ILP.input`, into individual `*.txt` files, one per question, that can be fed to `arithCplex`. The output is a set of files in `outdir` named `q000.txt, q001.txt, q002.txt,` and so on.

Usage: `splitQuestions.sh unsplitInputFile outdir`

Note that this is a very basic script that assumes a uniform format in the joint file, in terms of the number of lines (including empty lines) per question (default 9) and the number of lines to discard per question (currently 2). These may need to be adjusted manually if the joint input file format changes.

