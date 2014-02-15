load('../../sage_circuit_analysis.py')
trafo = SmallSignalLinearCircuit('trafo.cir', set_default_ic_to_zero=True)
symbolic_solution = trafo.solve_nodal_equations_symb()
numeric_solution = trafo.solve_nodal_equations_num()
