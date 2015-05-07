modes
=====

A multi-objective optimizer for the creation of peptide sequences mimicking
a set of reference peptides.

How does it work?
-----------------

`modes` allows you to generate a new joint peptide from a set of domains that
should be included in the new peptide and a list of reference peptides which are
ultimately used to create the joining peptide fragments. Thus, given an input set
of reference sequences `S` consiting of positive and negative examples to be mimicked
`modes` and a list of short peptide domains `D`=d_i, modes generates a set of fragments f_i
yielding the joined peptide `f_1-d_1-f_2-d_2-...-f_n-d_n-f_{n+1}` with the 
highest probability of belonging to `S`. The multi-objective capabilities allow you 
to optimize the fragments for a set of proteins types `S` simultaneously, where
the optimized quantity is the joint probability.

In order to achieve its goal `modes` combines classification using random forests 
on the reference sequences on a set of 26 physiochemical properties of the proteins 
with a gloabl optimization strategy using simulated annealing with energy landscape
pavement. Here `modes` runs all thos processes in parallel allowing large scale
design experiments.


Installation
------------

In order to compile `modes` you will need more-or-less recent C++ compiler
and Cmake (current versions of gcc or clang are perfect). For Linux cmake can
be installed with your package manager. `modes` is currently only tested on cmake.

In order to install git clone or unpack the downloaded zip of this repository.

A typical build could be run on Debian or Ubuntu like this for instance:

```bash
sudo apt-get install cmake build-essential git
cd ~/code
git clone https://github.com/cdiener/modes.git 
cd modes
mkdir bin && cd bin
cmake ..
make
```

This will install `modes` into you home directory under code/modes/bin. This will
be sufficient to run `modes`. From within that directory. To run it wherever you can
add it to you path with (add this line to your ~/.bashrc to make it permanent)

```bash
export PATH=PATH:~/code/modes/bin
``` 

If you additionally want to run the scripts reproducing
the analysis in the publication you will also need R and some of its packages.
To install R and all of its additional required packages use:

```bash
sudo apt-get install r-base r-base-dev
cd ~/code/modes/bin/scripts
./install_scripts
```

You can now run the scripts by hand or use `./run_all` to run all of them.

Usage
-----

## Design

In order to work `modes` requires the following things

1. A text (one sequence per line) or fasta file with the domains to be joined
2. One or more pairs of text or fasta files containing positive and negative 
   examples for the reference proteins. Those files should be named pos_NAME.txt
   and neg_NAME.txt or *.fasta. All pairs of those files should be located in a
   single directory.
3. For n domains you have to chose the n+1 maximum lengths of the joining 
   fragments k(1) to k(n+1).

If you have all those things you can call modes with
```bash
./modes domain_file 'k(1) ... k(n+1)' n_it ref_dir NAME1 NAME2 ...
``` 

Here, domain_file is the field containing you domains, the quoted k's are the 
maximum lengths of the joining fragments, n_it is the number of iterations for 
the optimizer, ref_dir is the path to the directory containingyour reference 
lists and NAMEk are the names of the (at least one) set of reference sequences. 
2.000 to 10.000 iterations are usually sufficient for convergence. 

This might sound complicated so here is some example where we try to design a peptide 
containing the nuclear localization signal (NLS) of the SV40 virus with the 
alpha-factor domain from yeast with in a way that the entire sequence mimicks a 
cell-penetrating peptide (CPP) with high efficiency. We also want no N-terminal 
odification of the alpha-domain and fragments of maximum length 4 elsewhere.All 
the files are provided in the example folder and we can run the entire work 
flow with:

```bash
./modes examples/alpha_NLS.txt '4 4 0' 2000 examles CPP efficiency
```

This will run 2000 iterations of optimization and show you the resulting
optimal fragments. The output will look something like this:

```
Classifying CPP on 27 variables over 1736 peptides................done.
Saved model and error estimates to examples/model_CPP.txt.
Training set classification error: 0.00633641
cv error estimate: 0.0924539 +- 0.00354885

Classifying efficiency on 27 variables over 374 peptides................done.
Saved model and error estimates to examples/model_efficiency.txt.
Training set classification error: 0.0160428
cv error estimate: 0.332194 +- 0.00998871

Needed 14.0503 s (+- 1e-09 s) for classification.

Optimizing on 8 threads.
Initial solution:
Iter: 0 (0 accepted) | Current: PKKKRKVWHWLQLKPGQPMY (E = 0.469) | best P(CPP): 0.531




Best sequences found: 
Links            P0      P1      P(all)  %α
NRKR RWKR _      0.915   0.86    0.787   0.6
NKRR RWKR _      0.915   0.86    0.787   0.6
NRKK RWRR _      0.915   0.86    0.787   0.6
RKKW MRRR _      0.935   0.84    0.785   0.6
RRWK MRKR _      0.935   0.84    0.785   0.6
KRWK MRRR _      0.935   0.84    0.785   0.6
RRWR MKKR _      0.935   0.84    0.785   0.6
RRKW MKRR _      0.935   0.84    0.785   0.6

Needed 2.87367 s (+- 1e-09 s) for optimization.
```

It shows you the errors of the initial classification training and finally the
best results with the respective fragments the probability for CPP (P0), the probability
for efficiency (P1) and the joint probability P(all) followed by a prediction for the
alpha-helical content of your peptide.

## Predict

In case you only want to run prediction on a set of sequences you can use the
`predict` program. Predict is called in the following way:

```bash
./predict data_dir 'C1 (C2) ...' data_file1 data_file2 ...
```

Here data_dir again is the directory with you positive and negative examples or
cached model files. C1 to Cn are the classes you want to predict. This is followed
by a list of data files in txt or fasta format. Predict will train a classifier
using the files in data_dir (or use a cached model if existing) and predict the average
P(C1,...,Cn) for each data file.

Contributors
------------

* Creator: Dr. Christian Diener (ch.diener[-at-]gmail.com)
