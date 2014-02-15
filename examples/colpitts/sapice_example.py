load('../../sage_circuit_analysis.py')
colpitts = SmallSignalLinearCircuit('colpitts_bjtmodel.cir', 'colpitts_bjtmodel.log')
symbolic_solution = colpitts.solve_nodal_equations_symb()
numeric_solution = colpitts.solve_nodal_equations_num()
