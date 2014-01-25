SAPICE
======

Symbolic analysis of electronic circuits in SAGE (sagemath) starting from a SPICE (ngspice) netlist.

Please see the documentation provided with the project and the documentation generated with doxygen (see corresponding directories).

This software will read a simple ngspice circuit netlist and compute the nodal equations for the corresponding linearized small circuit both symbolically and numerically.

The nodal equations can then be solved with the sage computer algebra system and it is possible to compute impedances and two-port network parameters both symbolically and numerically. Functions for automatic simplification of large expressions are provided.

Poles and zeroes can be computed.

The symbolic solution in this case will only be possible if the degree of the polynomial given by the denominator of the network function (in the Laplace domain) will be lesser than or equal to four.

Otherwise (symbolic) approximations are possible and of course a numerical solution can be computed.

Notice: the simple documentation file is at the moment out-of-date and needs to be rewritten. Please look at the code itself directly !

License: GNU GPL with absolutely no warranty !
