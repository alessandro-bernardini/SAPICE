load('../../sage_circuit_analysis.py')
rlc = SmallSignalLinearCircuit('RLC.cir', set_default_ic_to_zero=False)
symbolic_solution = rlc.solve_nodal_equations_symb()
numeric_solution = rlc.solve_nodal_equations_num()
