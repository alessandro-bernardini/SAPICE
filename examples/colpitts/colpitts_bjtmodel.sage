reset()
load("../../sage_circuit_analysis.py")
colpitts = SmallSignalLinearCircuit("colpitts_bjtmodel.cir", "colpitts_bjtmodel.log")
impedance5_0 = colpitts.impedance(('5','0'))
impedance5_0_num = impedance5_0.substitute(colpitts.dict_default_vals()[0])
impedance5_0_num_DC = impedance5_0_num.substitute({s:0})
sscolpitts_string = colpitts.export_lin_circuit_string(shorts={'3':'0', '1':'0', '4':'5'}, remove=[ 'LRFC' ])
sscolpitts = SmallSignalLinearCircuit(circuit_netlist=sscolpitts_string)
ssimpedance5_0 = sscolpitts.impedance(('5','0'))
dct = colpitts.dict_default_vals()[0]
simple_ssimpedance5_0 = simplify_rational_func(ssimpedance5_0, dct, 0.9)[0]
super_simple_ssimpedance5_0 = simplify_rational_func(ssimpedance5_0, dct, 1)[0]
