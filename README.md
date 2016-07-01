# RiboLogic

RiboLogic is a Python package for designing RNAs to adopt specified secondary structures in multiple conditions.

See `example_run.sh` for example usage.  Be sure to change `YOUR_CODE_DIR` to the path to your code.

# Dependencies

## Python

Python2.7+, Python packages: `numpy`, `matplotlib`, `requests`, `networkx`.

## NUPACK

Add the following to your bashrc
```shell
export NUPACKHOME=/path/to/RiboLogic/resources/nupack
```

# Input file format

Input files contain the following types of entries.

## Inputs

Inputs, including ligands and RNAs, are specified using the input symbol `<`.
For RNA inputs, specify the sequence on the subsequent line.

```
<myRNA
AUGCAUGC
```

For ligand inputs, specify the Kd of the ligand (in the same units as
concentration is specified subsequently) and the secondary structure
constraints of the aptamer

```
<myligand
1.5e-6
......................(((xxxxxx(((........................)))xxxxx))).....
```

## Sequence constraints

Sequence constraints are specified using the input symbol `-`.  Characters
must be `A`, `U`, `G`, `C` or `N`.

```
-sequence
NNNNNNNNNNAUGCNNNNNNNNNNN
```

Fixed subsequences that can be located anywhere in the full sequence can be
specified as follows.

```
-variable
AAAA
```

For each variable element, you must add two additional lines to the secondary
structure constraints specifying any constraints on the variable element.

In addition, you can specify subsequences that should not occur using `x`.

```
xGGGG
```

## Secondary structure constraints

Secondary structure constraints, indicated by `>`, may be specified for various conditions.  They
can have one of the following types: single (no inputs), oligos (RNA inputs),
and  aptamer (ligand inputs).  The constraints are specified by one line with
the desired secondary structure and a second line that specifies the bases
that should be constrained (`x` = constrained, `o` = unconstrained, `p` = paired (to any other base), `n` =
anything that is *not* what is given)

```
>single
.....(((((.........))))).
oooooxxxxxxxxxxxxxxxxxxxo
```

For oligos type, the second line should specify the RNAs and their concentrations (in M).

```
>oligos
myRNA 200e-6; myRNA2 100e-6
.....(((((.......))))))).
ooooooonnnooooooonnnnnooo
```

For aptamer type, the second line should specify the ligand and its
concentration.

```
>aptamer
myligand 200e-6
.....(((((.......))))))).
ooooooonnnooooooonnnnnooo
```

For bases that need to be paired but not to a specific base, use `|`.  You can supply an optional threshold parameter on the line following the constraints.  In the following example, the design must have at least 5 paired bases in the indicated region.

```
>single
.....||||||||............
oooooppppppppoooooooooooo
5
```
