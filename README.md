# AbaqusUsub.cpp

A framework to write Abaqus User (Material) Subroutine in C++.
For convenience, the FORTRAN vectors and matrices are reinterpreted into C++ vectors and matrices provided by the [LMT++](https://github.com/hleclerc/LMTpp) library.

### Installation

```
> git clone https://github.com/mgenet/LMTpp.git
> git clone https://github.com/mgenet/AbaqusUsub.cpp.git
```

### Usage

Examples are given for linear elasticity.
To compile them, simply run:

```
> cd AbaqusUsub.cpp
> make
```

### Example

An example is provided of a simple cube in tension.
To execute it, simply run:

```
> cd JOB
> make
```
