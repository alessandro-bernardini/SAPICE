load('../../sage_circuit_analysis.py')
lo_noise = SmallSignalLinearCircuit('lo_noise_osc.cir', 'lo_noise_osc.log', set_default_ic_to_zero=True)
symbolic_solution = lo_noise.solve_nodal_equations_symb()
numeric_solution = lo_noise.solve_nodal_equations_num()
