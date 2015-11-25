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

Basic usage: `arithCplex inputfile` where `inputfile` is a plain text file containing parameters to build an arithmetic expression model, such as the following in `example.txt`:

```
quantities : 8 3 x
types : "calorie" "bar" "calorie"
operators : + - * / =
n : 5
answer : 24.0
```

Note that the syntax of this input file is rather brittle. It uses white-space as delimiters and, in particular, expects a white-space on both sides of the `:` character.

The output includes a number of CPLEX related information as well as the number of solutions found. For each solution, it prints information in a format similar to the following excerpt from the output of `./arithCplex example.txt`:

```
parameters: n=5 l=2 k=1 p=5 q=2 m=8
TOTAL 25 solutions found
SOLN: CORRECT | POS/NEG | INT/FRA | OBJ-SCORE | TRUE-ANS | ANS | INFIX | POSTFIX | TYPED-POSTFIX
EXPR: 1 | POS | INT | 0 | 24 | 24 | (8*3)=x | 8 3 * x = | 8:calorie 3:bar *:calorie x:calorie =:calorie
EXPR: 1 | POS | INT | 0 | 24 | 24 | x=(8*3) | x 8 3 * = | x:calorie 8:calorie 3:bar *:calorie =:calorie
EXPR: 1 | POS | INT | 2 | 24 | 24 | 8=(x/3) | 8 x 3 / = | 8:calorie x:calorie 3:bar /:calorie =:calorie
EXPR: 1 | POS | INT | 4 | 24 | 24 | (x/3)=8 | x 3 / 8 = | x:calorie 3:bar /:calorie 8:calorie =:calorie
EXPR: 0 | POS | INT | 10 | 24 | 5 | (8-3)=x | 8 3 - x = | 8:calorie 3:calorie -:calorie x:calorie =:calorie
EXPR: 0 | POS | INT | 10 | 24 | 11 | (8+3)=x | 8 3 + x = | 8:calorie 3:calorie +:calorie x:calorie =:calorie
EXPR: 0 | POS | INT | 10 | 24 | 11 | x=(8+3) | x 8 3 + = | x:calorie 8:calorie 3:calorie +:calorie =:calorie
EXPR: 0 | POS | INT | 10 | 24 | 5 | x=(8-3) | x 8 3 - = | x:calorie 8:calorie 3:calorie -:calorie =:calorie
NET 8 unique, non-negative, integer-valued solutions found out of 25 total solutions
```

Convenience Utilities
---------------------

The `scripts` folder contains a couple of non-essential utilities.

**Running the CPLEX model on all `*.txt` files in a folder**

For experimentation, the included `runall.sh` script in the `scripts` folder may be used to run the CPLEX model on all `*.txt` files included in a specified folder. The output for each arithmetic problem file `file.txt` is then captured in `file.txt.out` in the same folder as `file.txt`.

Usage: `./runall.sh datafolder`


**Splitting a joint input file for a dataset into multiple `*.txt` files for `arithCplex`**

The `splitQuestions.sh` script in the `scripts` folder is a utility that may be used to split a joint input file, such as `tacl2015.ILP.input`, into individual `*.txt` files, one per question, that can be fed to `arithCplex`. The output is a set of files in `outdir` named `q000.txt, q001.txt, q002.txt,` and so on.

Usage: `splitQuestions.sh unsplitInputFile outdir`

Note that this is a very basic script that assumes a uniform format in the joint file, in terms of the number of lines (including empty lines) per question and the number of bottom and top lines to discard per question. These parameters may need to be adjusted manually if the joint input file format changes.

