#!/usr/bin/python
# -*- coding: utf-8 -*-

## @package sage_circuit_analysis
# \mainpage
# Author: Alessandro Bernardini
#
#     alessandro.bernardini@tum.de
#
# https://github.com/alessandro-bernardini/SAPICE
#
# http://alessandro-bernardini.github.io/SAPICE/
#
# License: GNU GPL
#
# Disclaimer: THERE IS NO WARRANTY FOR THE PROGRAM (SAPICE and all its provided 
# components), TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE
# STATED IN WRITING THE COPYRIGHT HOLDERS (Alessandro Bernardini) AND/OR OTHER 
# PARTIES PROVIDE THE PROGRAM AS IS WITHOUT WARRANTY OF ANY KIND, EITHER 
# EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS 
# TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM 
# PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR 
# CORRECTION.
#
# IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY 
# COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM 
# AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, 
# SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR 
# INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR 
# DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR 
# A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH 
# HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
#
# THIS PROGRAM WAS NOT WELL TESTED !
#
# Project home:
# https://github.com/alessandro-bernardini/SAPICE
#
# Requires: ngspice revision 24; sage version 5.6. Other version of ngspice or 
# sage should work well too.
#
# http://www.sagemath.org/
#
# http://ngspice.sourceforge.net/download.html
# 
# See documentation and license
#
#

"""@package docstring
Version: 0.3

Author: Alessandro Bernardini

alessandro.bernardini.tum@gmail.com    	
https://github.com/alessandro-bernardini/SAPICE
http://alessandro-bernardini.github.io/SAPICE/

License: GNU GPL.
http://www.gnu.org/licenses/gpl.html

Disclaimer: THERE IS NO WARRANTY FOR THE PROGRAM (SAPICE and all its provided 
components), TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE
STATED IN WRITING THE COPYRIGHT HOLDERS (Alessandro Bernardini) AND/OR OTHER 
PARTIES PROVIDE THE PROGRAM AS IS WITHOUT WARRANTY OF ANY KIND, EITHER 
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS 
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM 
PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR 
CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY 
COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM 
AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, 
SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR 
INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR 
DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR 
A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH 
HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

THIS PROGRAM WAS NOT WELL TESTED !

Project home:
https://github.com/alessandro-bernardini/SAPICE

Requires: ngspice revision 24; sage version 5.6. Other version of ngspice or 
sage should work well too.

http://www.sagemath.org/

http://ngspice.sourceforge.net/download.html

See documentation and license

"""

import sage.all as sage
import networkx as nx
import re
#import sys
import copy
#import string

RESISTOR_EXPR = re.compile('^(\s)*(R|r)\w+')
CAPACITOR_EXPR = re.compile('^(\s)*(C|c)\w+')
INDUCTOR_EXPR = re.compile('^(\s)*(L|l)\w+')
V_EXPR = re.compile('^(\s)*(V|v)\w+')
I_EXPR = re.compile('^(\s)*(I|i)\w+')
VCCS_EXPR = re.compile('^(\s)*(G|g)\w+')
BJT_EXPR = re.compile('^(\s)*(Q|q)\w+')
KIND_EXPR = re.compile('^(\s)*(K|k)\w+')
#DATA_FIELD = re.compile('\s+\w+\=-?\w*\.?\w+|\s+-?\w*\.?\w+|\s+\{\w+\}')
DATA_FIELD = re.compile('\s+\w+\=[-\+]?\w*\.?\w+[-\+]?\d*|\s+[-\+]?\w*\.?\w+[-\+]?\d*|\s+\{\w+\}')
TEMP_EXPR = re.compile('^(\s)*(\.temp)\s+\w+', re.IGNORECASE)


class SmallSignalLinearCircuit:
    """PROJECT HOME:

    	https://github.com/alessandro-bernardini/SAPICE
    	http://alessandro-bernardini.github.io/SAPICE/
    
    This class implements nodal analysis for linear electrical circuits.
    If a nonlinear circuit is given then the linearized small signal circuit is
    considered.
    This class is designed to be used WITHIN the sage computer algebra system

    	http://www.sagemath.org/

    and in conjunction with ngspice (open source version of spice):

    	http://ngspice.sourceforge.net/

    So you have to install sage and ngspice first (see relative documentation).

    Then run	

    	sage

    from the shell
    then in the sage prompt import the present module

    	import sage_circuit_analysis

    and the create a SmallSignalLinearCircuit object

    	circuit = sage_circuit_analysis.SmallSignalLinearCircuit('circuitfile.cir')

    with circuitfile.cir a spice netllist (see additional information below
    and in the help of the constructor).

    You can read-in a ngspice netlist and typically you will provide the operating
    point data and other information via an ngspice batch output file
    that you have to generate with the command (from the shell):
    
    	ngspice -b circuitfile.cir -o batchoutput.log

    where circuitfile.cir is the spice netlist
    and batchoutput.log is the log file containing the output data of interest.
    Both files circuitfile.cir and batchoutput.log must be passed
    to a SmallSignalLinearCircuit object when invoking the constructor.

    The nodal equations can be solved both numerically (if numerical values
    are provided) and symbolically.

    You can compute impedances and two port network parameters and do
    the computation of poles and zeroes if a symbolic closed form solutions
    exists.

    It is possible to compute symbolic approximations to the exact solution
    where only dominant terms are considered.

    NOTE: numerical solutions usually means a polynom or a rational function
    dependent on the complex variable s and on independent sources.

    Usually all quantities are complex.
    """

    # when implementing modified nodal analysis a set of MNA_equations will 
    # be needed togheter with additional MNA_unknowns, that usually will 
    # represent currents in branches (edges)

    def __init__(
        self,
        filename=None, # a filename containing a circuit netlist
        spice_batch_output_file=None,
        circuit_netlist=None,#string describing a LINEAR ngspice circuit netlist 
        check_operating_region=True,
        set_default_ic_to_zero=True,
        ignore_all_ic=False,
        ):
        """filename='circuitfile.cir'
	where circuitfile.cir is a ngspice (spice) netlist.

	spice_batch_output_file='batchoutput.log'
	a file needed when operating point data must be considered
	You have to run from the shell the command
	ngspice -b circuitfile.cir -o batchoutput.log
	fog generating the batchoutput.log file

	in alternative
	circuit_netlist=CIRCUIT_STRING
	can be used where CIRCUIT_STRING is a string containing
	the circuit netlist (in place of circuitfile.cir)

        with check_operating_region=True a linearized model in dependence of the 
        operating region will be choosen for each semiconductor devices. 
        Otherwise it will be assumed that a default model is valid. (for 
        BJT transistors the active region small signal model will be this 
        default model)
    
        ignore_all_IC=True will ignore all initial conditions for both symbolic
        and numeric computations
    
        set_default_ic_to_zero=True will set initial conditions to zero for
        capacitors or inductors where no initial condition is specified
        """

        if filename != None:
            if circuit_netlist != None:
                raise Exception('You cannot give both a circuit file and a '
                                'linear circuit netlist in string form when '
                                'constructing a SmallSignaLinearCircuit.'
                                )
            else:
                try:
                    main_file = open(filename, 'r')
                except:
                    raise IOError(filename + ' file does not exist')
                self.circuit_file = main_file.readlines()
                self.original_circuit_file = self.circuit_file
                if spice_batch_output_file != None:
                    try:
                        bof = open(spice_batch_output_file, 'r')
                    except:
                        raise IOError(spice_batch_output_file
                                + ' file does not exist')
                    self.spice_batch_output_file_original = bof.readlines()
                else:
                    self.spice_batch_output_file_original = []
        else:
            if circuit_netlist == None:
                raise Exception('Provide a circuit file or a linear circuit '
                                'netlist string when constructing a '
                                'SmallSignaLinearCircuit object.'
                                )
            else:

                # circuit_netlist != None

                self.circuit_file = circuit_netlist.splitlines(True)
                self.original_circuit_file = \
                    circuit_netlist.splitlines(True)

    # initialize to empty entities and set default temperature.

        self.circuit_graph = nx.MultiDiGraph()
        self.circuit_graph_substitutions = None  # will be generated later
        self.sources = []
        self.sources_names = []
        self.sources_component_id = []
        self.circuit_variables = []
        self.circuit_parameters = []
        self.coupled_inductors_matrix = {}
        self.couplings = {}
        self.coupled_inductors_data = {}
        self.default_substitutions = {}  # stores the default substitutions with string values
        self.default_substitutions_values = {sage.var('V_0'): 0}  # stores the values for the default substitutions; initialized with V_0 = 0
        self.nodal_equations = {}  # currents flowing out of a node are considered to be positive.
        # self.nodal_equations is a dict. The keys are the node names
        # the values are the KCL equations for the corresponding node
        # currents flowing OUT of a given node are considered POSITIVE
        # in the corresponding KCL equation (the value in the dict)
        # for the considered node (the key in the dict)
        self.nodal_equations_substitutions = []
        self.additional_equations = []
        self.additional_equations_explicit = {}
        self.additional_equations_exp_sub = {}
        self.initial_conditions = []
        self.thermal_voltages = {}
        self.nodal_voltages = []
        self.I_cpld_ind_solutions = None
        self.V_cpld_ind = {}
        self.I_cpld_ind = {}
        self.M_cpld = {}
	self.BJT_unknowns = []
        self.operating_regions = {}  # remember the operating regions for the semiconductor devices
        self.temp = 300.15  # default value. Can be overridden if a .temp line is present
        pos = 0
	self.spice_batch_output_file = self._preprocessngspicelog(self.spice_batch_output_file_original)

# first for loop over self.circuit_file lines
# in this part preliminary preparatives are taken for a later nodal analysis.


#------------------------------------------------------------------------ 
        #first for loop over lines in circuit file for handling BJT, TEMP
        for line in self.circuit_file:
            pos = pos + 1

            if TEMP_EXPR.match(line):
                data_temp = DATA_FIELD.findall(line)  # will include .temp at index 0
                try:
                    temperature = data_temp[0].strip()
                    self.temp = extract_value(temperature) + 273.15  
                    # sets the temperature of the circuit in kelvin
                except:
                    self.temp = 300.15
                    print 'Error in setting the temperature. Default value of 27 C = 300.15 K is assumed'

            if V_EXPR.match(line):
                raise Exception('Voltage sources are for now not supported. '
                                'Replace them with an electrically equivalent'
                                ' current source (source transformation)'
                                )

            if BJT_EXPR.match(line):
                bjt_lineptr = BJT_EXPR.match(line)

        # bjt_line_data stores the data of the bjt line

                bjt_line_data = line[bjt_lineptr.end():]
                bjt_id = line[:bjt_lineptr.end()].strip()
                operating_region = 'unknown'

            # if check_operating_region is not set then it will be assumed that the bjt bjt_id operates in its active region

                data_bjt = DATA_FIELD.findall(bjt_line_data)
                data_bjt_strip = []

                # data_bjt_strip describes the (stripped) data fields of the bjt line line

                for key_neq in data_bjt:
                    data_bjt_strip = data_bjt_strip + [key_neq.strip()]

                # data_bjt_strip contains the stripped strings representing the data fields relative to the line line describing a bjt in the spice file
                # the first 3 data fields describe the nodes to which the bjt is connected

                for key_neq in range(3):
                    if not sage.var('V_' + data_bjt_strip[key_neq]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_bjt_strip[key_neq])]
                        if data_bjt_strip[key_neq] != '0' \
                            and data_bjt_strip[key_neq] != 'gnd' \
                            and data_bjt_strip[key_neq] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_'
                                   + data_bjt_strip[key_neq])]

                try:
                    model = data_bjt_strip[3]
                except IndexError:
                    raise Exception('Provide transistor model for transistor connected to nodes '
                                     + data_bjt_strip[0] + ' '
                                    + data_bjt_strip[1] + ' '
                                    + data_bjt_strip[2]
                                    + ' (respectively C B E nodes) ')
                bjt_model_desc_reached = False
                model_present_bjt_reached = False
                skipnextline = False
                for key_neq in self.spice_batch_output_file:
                    if skipnextline == False:
                        if re.match('^\sBJT\smodels', key_neq) != None:
                            bjt_model_desc_reached = True
                            skipnextline = False
                        elif re.match('^\s+model\s+' + model, key_neq,
                                flags=re.IGNORECASE) != None:

                            model_present_bjt_reached = True
                            skipnextline = True
                        elif bjt_model_desc_reached \
                            and model_present_bjt_reached:
                            matchline = \
                                re.match('^\s+(?P<field>[a-zA-Z]+)\s+(?P<value>(((\+|-)?\d*\.?\d+(e(\+|-)?\d+)?)|(npn)|(pnp)))'
                                    , key_neq, re.IGNORECASE)
                            if matchline != None:
                                if matchline.group('field') != 'type':

                                    self.default_substitutions[sage.var(bjt_id
        + '_' + matchline.group('field').upper())] = \
    matchline.group('value')
                                else:
                                    if not matchline.group('value') \
    == 'npn':
                                        print 'Only npn bjt are supported '
                                        print matchline.group('field') \
    + ' ' + matchline.group('value') + ' is not supported'
                                        print matchline.group(0)
                                        raise Exception('Device is not supported'
        )
                            skipnextline = False
                            if re.match('^\s*\n', key_neq) != None:
                                bjt_model_desc_reached = False
                                model_present_bjt_reached = False
                    else:
                        skipnextline = False
                bjt_reached = False
                present_bjt_reached = False
                skipnextline = False
                for key_neq in self.spice_batch_output_file:
                    if skipnextline == False:
                        if re.match('^\sBJT:\sBipolar', key_neq) \
                            != None:
                            bjt_reached = True
                            skipnextline = False
                        elif re.match('^\s+device\s+' + bjt_id,
                                key_neq, flags=re.IGNORECASE) != None:

                            present_bjt_reached = True
                            skipnextline = False
                        elif bjt_reached and present_bjt_reached:
                            matchline = \
                                re.match('^\s+(?P<field>[a-zA-Z]+)\s+(?P<value>(\+|-)?\d*\.?\d+(e(\+|-)?\d+)?)'
                                    , key_neq, re.IGNORECASE)
                            if matchline != None:
                                if matchline.group('field') in ['vbe',
                                        'vbc', 'ib', 'ic']:

                                    self.default_substitutions[sage.var(bjt_id
        + '_' + matchline.group('field').upper() + '_Q')] = \
    matchline.group('value')
                            skipnextline = False
                        elif re.match('^\s*\n', key_neq) != None:
                            bjt_reached = False
                            present_bjt_reached = False
                            skipnextline = False
                    else:
                        skipnextline = False
                try:
                    if extract_value(self.default_substitutions[sage.var(bjt_id
                            + '_VBC_Q')]) < 0 \
                        and extract_value(self.default_substitutions[sage.var(bjt_id
                            + '_VBE_Q')]) > 0:
                        operating_region = 'active'
                        self.operating_regions[bjt_id] = 'active'
                    else:
                        operating_region = 'unknown'
                except:
                    operating_region = 'unknown'

                if check_operating_region and operating_region \
                    == 'active' or not check_operating_region:
                    self.circuit_parameters = self.circuit_parameters \
                        + [
                        sage.var(bjt_id + '_BFop'),
                        sage.var(bjt_id + '_IB_Q'),
                        sage.var(bjt_id + '_VCE_Q'),
                        sage.var(bjt_id + '_VAF'),
                        sage.var(bjt_id + '_IC_Q'),
                        sage.var(bjt_id + '_VT'),
                        sage.var(bjt_id + '_TF'),
                        sage.var(bjt_id + '_CJC'),
                        sage.var(bjt_id + '_VBC_Q'),
                        sage.var(bjt_id + '_VJC'),
                        sage.var(bjt_id + '_MJC'),
                        sage.var(bjt_id + '_VBE_Q'),
                        ]

            # add equation relating vce with vbe and vbc. 
            # This equation is valid in any case, but we still introduce it 
            # here in the additional_equations....

                    self.additional_equations = \
                        self.additional_equations + [sage.var(bjt_id
                            + '_VCE_Q') == sage.var(bjt_id + '_VBE_Q')
                            - sage.var(bjt_id + '_VBC_Q')]

            # add equation describing the output resistance of the linearized 
            # circuit of the BJT

                    self.additional_equations = \
                        self.additional_equations + [sage.var('RO_'
                            + bjt_id) == (sage.var(bjt_id + '_VCE_Q')
                            + sage.var(bjt_id + '_VAF'))
                            / sage.var(bjt_id + '_IC_Q')]

            # add equation describing the transconductance Gm

                    self.additional_equations = \
                        self.additional_equations + [sage.var('GM_'
                            + bjt_id) == sage.var(bjt_id + '_IC_Q')
                            / sage.var(bjt_id + '_VT')]

            # add equation describing the capacitance between base and emitter 
            # in the linearized circuit

                    self.additional_equations = \
                        self.additional_equations + [sage.var('CBE_'
                            + bjt_id) == sage.var(bjt_id + '_TF')
                            * sage.var('GM_' + bjt_id)]

            # add equation describing the capacitance between base and 
            # collector in the linearized circuit

                    self.additional_equations = \
                        self.additional_equations + [sage.var('CBC_'
                            + bjt_id) == sage.var(bjt_id + '_CJC') * (1
                            - sage.var(bjt_id + '_VBC_Q')
                            / sage.var(bjt_id + '_VJC'))
                            ** -sage.var(bjt_id + '_MJC')]

            # add equation describing input resistance of the linearized circuit

                    self.additional_equations = \
                        self.additional_equations + [sage.var('RPI_'
                            + bjt_id) == sage.var(bjt_id + '_VT')
                            / sage.var(bjt_id + '_IB_Q')]

            # add equation describing beta of the transistor

                    self.additional_equations = \
                        self.additional_equations + [sage.var(bjt_id
                            + '_BFop') == sage.var(bjt_id + '_IC_Q')
                            / sage.var(bjt_id + '_IB_Q')]
                    self.additional_equations = \
                        self.additional_equations + [sage.var(bjt_id
                            + '_VT') == sage.var('GLOBAL_K_BOLTZMANN')
                            * sage.var('GLOBAL_TEMP')
                            / sage.var('GLOBAL_q')]
                    self.circuit_file = self.circuit_file[:pos - 1] \
                        + ['RO_' + bjt_id + '  ' + data_bjt_strip[0]
                           + ' ' + data_bjt_strip[2] + ' {RO_' + bjt_id
                           + '} \t\t\n'] + ['GM_' + bjt_id + '  '
                            + data_bjt_strip[0] + ' '
                            + data_bjt_strip[2] + ' '
                            + data_bjt_strip[1] + ' '
                            + data_bjt_strip[2] + ' {GM_' + bjt_id
                            + '}  \t\t\n'] + ['CBC_' + bjt_id + '  '
                            + data_bjt_strip[0] + ' '
                            + data_bjt_strip[1] + ' {CBC_' + bjt_id
                            + '} \t\t\n'] + ['CBE_' + bjt_id + '  '
                            + data_bjt_strip[1] + ' '
                            + data_bjt_strip[2] + ' {CBE_' + bjt_id
                            + '} \t\t\n'] + ['RPI_' + bjt_id + '  '
                            + data_bjt_strip[1] + ' '
                            + data_bjt_strip[2] + ' {RPI_' + bjt_id
                            + '} \t\t\n'] + self.circuit_file[pos:]
                    try:
                        '''add_eq_sol = \
                            sage.solve(self.additional_equations, [
                            sage.var('RO_' + bjt_id),
                            sage.var('GM_' + bjt_id),
                            sage.var('CBE_' + bjt_id),
                            sage.var('CBC_' + bjt_id),
                            sage.var('RPI_' + bjt_id),
                            sage.var(bjt_id + '_BFop'),
                            sage.var(bjt_id + '_VCE_Q'),
                            sage.var(bjt_id + '_VT'),
                            ], solution_dict=True)
			'''    
			self.BJT_unknowns += [sage.var('RO_' + bjt_id),\
                            sage.var('GM_' + bjt_id),\
                            sage.var('CBE_' + bjt_id),\
                            sage.var('CBC_' + bjt_id),\
                            sage.var('RPI_' + bjt_id),\
                            sage.var(bjt_id + '_BFop'),\
                            sage.var(bjt_id + '_VCE_Q'),\
                            sage.var(bjt_id + '_VT')]
			add_eq_sol = sage.solve(self.additional_equations, self.BJT_unknowns, solution_dict=True)

                        if len(add_eq_sol) == 1:
                            self.additional_equations_explicit.update(add_eq_sol[0])
                        else:
                            raise RuntimeError
                    except:
                        raise Exception('An error has occurred in solving the additional equations for '
                                 + bjt_id)
                else:

                    raise Exception('transistor ' + bjt_id
                                    + ' is not operating in a supported region'
                                    )

#-------------------------------------------------------
        #second for loop for handling COUPLED INDUCTORS:
        for line in self.circuit_file:
            if KIND_EXPR.match(line):
                kind_lineptr = KIND_EXPR.match(line)
                kind_line_data = line[kind_lineptr.end():]
                kind_id = line[:kind_lineptr.end()].strip()
                data_kind = DATA_FIELD.findall(kind_line_data)
                data_kind_strip = []
                for key_neq in data_kind:
                    data_kind_strip = data_kind_strip + [key_neq.strip()]
                try:
                    self.coupled_inductors_matrix[data_kind_strip[0]].update({data_kind_strip[1] : kind_id})
                except KeyError:
                    self.coupled_inductors_matrix[data_kind_strip[0]] = {data_kind_strip[1] : kind_id}
                try:
                    self.coupled_inductors_matrix[data_kind_strip[1]].update({data_kind_strip[0] : kind_id})
                except KeyError:
                    self.coupled_inductors_matrix[data_kind_strip[1]] = {data_kind_strip[0] : kind_id}

                try:
                    self.couplings[kind_id] = {'coupled_inductors': [None, None], 
                                                            'coupling_coefficient': data_kind_strip[2]}
                    self.circuit_variables += [ sage.var(kind_id) ]
                    self.default_substitutions[sage.var(kind_id)] = data_kind_strip[2]
                except IndexError:
                    self.couplings[kind_id] = {'coupled_inductors': [None, None], 
                                                            'coupling_coefficient': '{'+kind_id+'}'}
                    self.circuit_variables += sage.var(kind_id)

                for line2 in self.circuit_file:
                    #nested for loop for finding the inductors that are coupled via K
                    if INDUCTOR_EXPR.match(line2):
                        ind_lineptr = INDUCTOR_EXPR.match(line2)
                        ind_id = line2[:ind_lineptr.end()].strip()
                        self.coupled_inductors_data[ind_id] = \
                            get_inductor_data(line2,
                                              ignore_all_ic=ignore_all_ic,
                                              set_default_ic_to_zero=set_default_ic_to_zero)
                                                                                
                        if ind_id == data_kind_strip[0]:
                            self.couplings[kind_id]['coupled_inductors'][0] = line2       
                        elif ind_id == data_kind_strip[1]:
                            self.couplings[kind_id]['coupled_inductors'][1] = line2
        for ind_id in self.coupled_inductors_data.keys():
            if not ind_id in self.coupled_inductors_matrix.keys():
                self.coupled_inductors_data.pop(ind_id)
        #coupled_inductor_matrix and coupling factors constructed.
        
        for ind_id in self.coupled_inductors_data.keys():
            self.circuit_variables += self.coupled_inductors_data[ind_id]['circuit_variables']
            self.nodal_voltages += self.coupled_inductors_data[ind_id]['nodal_voltages']
            self.initial_conditions += self.coupled_inductors_data[ind_id]['initial_conditions']
            self.default_substitutions.update(self.coupled_inductors_data[ind_id]['default_substitutions'])
            try:
                self.circuit_graph.add_edge(self.coupled_inductors_data[ind_id]['node0'],
                                            self.coupled_inductors_data[ind_id]['node1'],
                                            type='L', id=self.coupled_inductors_data[ind_id]['id'],
                                            coupled_to={'coupled_inductors':self.coupled_inductors_matrix[ind_id],
                                                        'coupling_coefficient_data':self.couplings},
                                            value=self.coupled_inductors_data[ind_id]['value_in_line'])
            except:
                raise RuntimeError("An error has occurred")
        #beginning of computation for voltages on coupled inductors and for currents in nodal equations for coupled inductors
        self.I_cpld_ind = {}
        self.I_cpld_ind_list = []
        self.V_cpld_ind = {}
        self.M_cpld = {}
        for inductor_init in self.coupled_inductors_matrix.keys():
            self.V_cpld_ind[inductor_init] = 0
            self.I_cpld_ind[inductor_init] = None
        # initialization to None
        for cpld_inductor in self.coupled_inductors_matrix.keys():
            #compute voltages equations and then nodal equations
            cpld_inductor_data = self.coupled_inductors_data[cpld_inductor]
            cpld_inductor_node = {}
            cpld_inductor_node[0] = cpld_inductor_data['node0']
            cpld_inductor_node[1] = cpld_inductor_data['node1']
            try:
                K_coeff = 0
                for inductor in self.coupled_inductors_matrix[cpld_inductor].keys():
                    inductor_node0 = self.coupled_inductors_data[inductor]['node0']
                    inductor_node1 = self.coupled_inductors_data[inductor]['node1']
                    K_coeff = sage.var(self.coupled_inductors_matrix[cpld_inductor][inductor])
                    try:
                        self.M_cpld[cpld_inductor].update({inductor: K_coeff*sage.sqrt(sage.var(cpld_inductor)\
                                                                             *sage.var(inductor))})
                    except KeyError:
                        self.M_cpld[cpld_inductor] = {inductor: K_coeff*sage.sqrt(sage.var(cpld_inductor)\
                                                                             *sage.var(inductor))}
                    self.V_cpld_ind[cpld_inductor] = self.V_cpld_ind[cpld_inductor] \
                        + sage.var('s')*self.M_cpld[cpld_inductor][inductor]*sage.var('I_'+inductor)\
                        - self.M_cpld[cpld_inductor][inductor]*sage.var('I_initial_'+inductor +'_'\
                                                                       + inductor_node0 + '_'\
                                                                       + inductor_node1)
                    if self.I_cpld_ind[inductor] == None:
                        self.I_cpld_ind[inductor] = (sage.var('I_'+inductor), inductor_node0,
                                                     inductor_node1)
                        self.I_cpld_ind_list += [ sage.var('I_'+inductor) ]
                self.V_cpld_ind[cpld_inductor] = self.V_cpld_ind[cpld_inductor] \
                    + sage.var('s')*sage.var(self.coupled_inductors_data[cpld_inductor]['id'])*\
                    sage.var('I_'+cpld_inductor) - \
                    sage.var(self.coupled_inductors_data[cpld_inductor]['id'])*\
                    sage.var('I_initial_'+cpld_inductor +'_'+ cpld_inductor_node[0] +'_'+
                             cpld_inductor_node[1])
                if self.I_cpld_ind[cpld_inductor] == None:
                    self.I_cpld_ind[cpld_inductor] = (sage.var('I_'+cpld_inductor),
                                                      cpld_inductor_node[0],
                                                      cpld_inductor_node[1])
                    self.I_cpld_ind_list += [ sage.var('I_'+cpld_inductor) ]
            except:
                raise Exception("An error has occurred")

        temp_dict = {}
        for eqn in self.V_cpld_ind.keys():
            temp_dict[eqn] = self.V_cpld_ind[eqn] == sage.var('V_' + eqn)
        self.V_cpld_ind = temp_dict
        # voltage equations across inductors are constructed
        # the signs follow the dot convention: see ngspice documentation on coupled 
        # mutual inductors
        #voltage equations for coupled inductors are stated: now thay has to be solved
        self.I_cpld_ind_solutions = sage.solve(self.V_cpld_ind.values(), self.I_cpld_ind_list)
        if len(self.I_cpld_ind_list) == 1:
            raise WrongUse("If coupled inductors are used there must be at least two mutually coupled inductors")
        cpld_voltage_subst = {}
        for inductor in self.coupled_inductors_matrix.keys():
            cpld_voltage_subst[sage.var('V_'+inductor)] = sage.var('V_' + self.coupled_inductors_data[inductor]['node0']) \
                - sage.var('V_' + self.coupled_inductors_data[inductor]['node1'])
        simpl_temp = []
        for inductor_current in self.I_cpld_ind_solutions[0]:
            simpl_temp += [inductor_current.full_simplify()]
        self.I_cpld_ind_solutions[0] = simpl_temp
        for inductor in self.coupled_inductors_matrix.keys():
            for inductor_current in self.I_cpld_ind_solutions[0]:
                if (inductor_current.lhs() == sage.var('I_'+inductor)).test_relation():
                    if self.coupled_inductors_data[inductor]['node0'] not in ['0', 'gnd', 'GND']:
                        try:
                            self.nodal_equations[self.coupled_inductors_data[inductor]['node0']] = \
                                self.nodal_equations[self.coupled_inductors_data[inductor]['node0']] \
                                + inductor_current.rhs().substitute(cpld_voltage_subst)
                        except KeyError:
                             self.nodal_equations[self.coupled_inductors_data[inductor]['node0']] = \
                                inductor_current.rhs().substitute(cpld_voltage_subst)
                    if self.coupled_inductors_data[inductor]['node1'] not in ['0', 'gnd', 'GND']:
                        try:
                            self.nodal_equations[self.coupled_inductors_data[inductor]['node1']] = \
                                self.nodal_equations[self.coupled_inductors_data[inductor]['node1']] \
                                - inductor_current.rhs().substitute(cpld_voltage_subst)
                        except KeyError:
                             self.nodal_equations[self.coupled_inductors_data[inductor]['node1']] = \
                                - inductor_current.rhs().substitute(cpld_voltage_subst)
                     
        #nodal equations for coupled inductors branches are stated.
        #coupled inductors handled
                            

#-------------------------------------------------------
        # third  for loop  over self.circuit_file lines
        # proper nodal analysis done here.
        # nodal equations constructed here
        for line in self.circuit_file:
                           
            if I_EXPR.match(line):
                i_lineptr = I_EXPR.match(line)

            # line describes an independent current source in the spice circuit 
            # netlist file

                i_line_data = line[i_lineptr.end():]
                i_id = line[:i_lineptr.end()].strip()
                data_i = DATA_FIELD.findall(i_line_data)
                data_i_strip = []

        # data_i_strip describes the (stripped) data fields of the current 
        # source line line.

                for key_neq in data_i:
                    data_i_strip = data_i_strip + [key_neq.strip()]

        # data_i_strip contains the stripped strings representing the data 
        # fields relative to the line line describing a current source in the 
        # spice file
        
        # the first two data fields describes the nodes n+ and n- to which the 
        # current source is connected

                for key_neq in range(2):
                    if not sage.var('V_' + data_i_strip[key_neq]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_i_strip[key_neq])]
                        if data_i_strip[key_neq] != '0' \
                            and data_i_strip[key_neq] != 'gnd' \
                            and data_i_strip[key_neq] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_'
                                   + data_i_strip[key_neq])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(i_id)]

        # only the small signal offset will be considered as a current source. 
        # For DC currents (or DC components) a current equal to zero is 
        # considered. In case of only a DC current no modification to the nodal
        # equations is needed.

                self.circuit_graph.add_edge(
                    data_i_strip[0],
                    data_i_strip[1],
                    type='I',
                    id=i_id,
                    value=data_i_strip[2:],
                    ngspice_line_data=data_i_strip[0:2],
                    )
                try:
                    if data_i_strip[0] != '0' and data_i_strip[0] \
                        != 'gnd' and data_i_strip[0] != 'GND':
                        self.nodal_equations[data_i_strip[0]] = \
                            self.nodal_equations[data_i_strip[0]] \
                            + sage.function(i_id + '_small_signal_'
                                + data_i_strip[0] + '_'
                                + data_i_strip[1], sage.var('s'))
                except KeyError:
                    if data_i_strip[0] != '0' and data_i_strip[0] \
                        != 'gnd' and data_i_strip[0] != 'GND':
                        self.nodal_equations[data_i_strip[0]] = \
                            +sage.function(i_id + '_small_signal_'
                                + data_i_strip[0] + '_'
                                + data_i_strip[1], sage.var('s'))
                try:
                    if data_i_strip[1] != '0' and data_i_strip[1] \
                        != 'gnd' and data_i_strip[1] != 'GND':
                        self.nodal_equations[data_i_strip[1]] = \
                            self.nodal_equations[data_i_strip[1]] \
                            - sage.function(i_id + '_small_signal_'
                                + data_i_strip[0] + '_'
                                + data_i_strip[1], sage.var('s'))
                except KeyError:
                    if data_i_strip[1] != '0' and data_i_strip[1] \
                        != 'gnd' and data_i_strip[1] != 'GND':
                        self.nodal_equations[data_i_strip[1]] = \
                            -sage.function(i_id + '_small_signal_'
                                + data_i_strip[0] + '_'
                                + data_i_strip[1], sage.var('s'))

        # note: ngspice syntax is IYYYYYYY N+ N-
        # a current source of positive value forces current to flow out of n+ 
        # node, through the current source and into the n- node
        # n+ is associated to dara_i_strip[0] and n- to data_i_strip[1]

                self.sources += [sage.function(i_id + '_small_signal_'
                                 + data_i_strip[0] + '_'
                                 + data_i_strip[1], sage.var('s'))]
                self.sources_names += [i_id + '_small_signal_'
                        + data_i_strip[0] + '_' + data_i_strip[1]]
                self.sources_component_id += [i_id]

        # i_id + '_small_signal_' + data_i_strip[0] + '_' + data_i_strip[1] is 
        # the current imposed by source i_id flowing from data_i_strip[0] 
        # to data_i_strip[1] through the current source i_id itself.

            if RESISTOR_EXPR.match(line):
                res_lineptr = RESISTOR_EXPR.match(line)

                # line describes a resistor in the spice circuit netlist file

                res_line_data = line[res_lineptr.end():]
                resistor_id = line[:res_lineptr.end()].strip()
                data_res = DATA_FIELD.findall(res_line_data)
                data_res_strip = []

                # data_res_strip describes the (stripped) data fileds of 
                # the resistor line line.

                for key_neq in data_res:
                    data_res_strip = data_res_strip + [key_neq.strip()]

                # data_res_strip contains the stripped strings representing the
                # data fileds relative to the line line describing a resistor 
                # in the spice file
                
                # the first two data fields describes the nodes to which 
                #the resistor is connected

                for key_neq in range(2):
                    if not sage.var('V_' + data_res_strip[key_neq]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_res_strip[key_neq])]
                        if data_res_strip[key_neq] != '0' \
                            and data_res_strip[key_neq] != 'gnd' \
                            and data_res_strip[key_neq] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_'
                                   + data_res_strip[key_neq])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(resistor_id)]
                try:
                    self.default_substitutions[sage.var(resistor_id)] = \
                        data_res_strip[2]  
                    self.circuit_graph.add_edge(data_res_strip[0],
                            data_res_strip[1], type='R',
                            id=resistor_id, value=data_res_strip[2:]) 
                except IndexError:
                    self.circuit_graph.add_edge(data_res_strip[0],
                            data_res_strip[1], type='R', id=resistor_id)
                try:
                    if data_res_strip[0] != '0' and data_res_strip[0] \
                        != 'gnd' and data_res_strip[0] != 'GND':
                        self.nodal_equations[data_res_strip[0]] = \
                            self.nodal_equations[data_res_strip[0]] \
                            + (sage.var('V_' + data_res_strip[0])
                               - sage.var('V_' + data_res_strip[1])) \
                            / sage.var(resistor_id)
                except KeyError:
                    if data_res_strip[0] != '0' and data_res_strip[0] \
                        != 'gnd' and data_res_strip[0] != 'GND':
                        self.nodal_equations[data_res_strip[0]] = \
                            (sage.var('V_' + data_res_strip[0])
                             - sage.var('V_' + data_res_strip[1])) \
                            / sage.var(resistor_id)
                try:
                    if data_res_strip[1] != '0' and data_res_strip[1] \
                        != 'gnd' and data_res_strip[1] != 'GND':
                        self.nodal_equations[data_res_strip[1]] = \
                            self.nodal_equations[data_res_strip[1]] \
                            + (sage.var('V_' + data_res_strip[1])
                               - sage.var('V_' + data_res_strip[0])) \
                            / sage.var(resistor_id)
                except KeyError:
                    if data_res_strip[1] != '0' and data_res_strip[1] \
                        != 'gnd' and data_res_strip[1] != 'GND':
                        self.nodal_equations[data_res_strip[1]] = \
                            (sage.var('V_' + data_res_strip[1])
                             - sage.var('V_' + data_res_strip[0])) \
                            / sage.var(resistor_id)
            elif CAPACITOR_EXPR.match(line):
                cap_lineptr = CAPACITOR_EXPR.match(line)
                cap_line_data = line[cap_lineptr.end():]
                capacitor_id = line[:cap_lineptr.end()].strip()
                data_cap = DATA_FIELD.findall(cap_line_data)
                data_cap_strip = []

                # data_cap_strip describes the (stripped) data fileds of 
                # the capacitor line line.

                for key_neq in data_cap:
                    data_cap_strip = data_cap_strip + [key_neq.strip()]

                # data_cap_strip contains the stripped strings representing 
                # the data fileds relative to the line line describing 
                # a capacitor in the spice file
                # the first two data fields describes the nodes to which 
                # the capacitor is connected

                for key_neq in range(2):
                    if not sage.var('V_' + data_cap_strip[key_neq]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_cap_strip[key_neq])]
                        if data_cap_strip[key_neq] != '0' \
                            and data_cap_strip[key_neq] != 'gnd' \
                            and data_cap_strip[key_neq] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_'
                                   + data_cap_strip[key_neq])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var('V_initial_' + capacitor_id + '_'
                       + data_cap_strip[0] + '_' + data_cap_strip[1])]
                self.initial_conditions += [sage.var('V_initial_'
                        + capacitor_id + '_' + data_cap_strip[0] + '_'
                        + data_cap_strip[1])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(capacitor_id)]
                try:
                    self.default_substitutions[sage.var(capacitor_id)] = \
                        data_cap_strip[2]  
                    self.circuit_graph.add_edge(data_cap_strip[0],
                            data_cap_strip[1], type='C',
                            id=capacitor_id, value=data_cap_strip[2:])  
                except IndexError:
                    self.circuit_graph.add_edge(data_cap_strip[0],
                            data_cap_strip[1], type='C',
                            id=capacitor_id)
                try:
                    if data_cap_strip[0] != '0' and data_cap_strip[0] \
                        != 'gnd' and data_cap_strip[0] != 'GND':
                        self.nodal_equations[data_cap_strip[0]] = \
                            self.nodal_equations[data_cap_strip[0]] \
                            + (sage.var('V_' + data_cap_strip[0])
                               - sage.var('V_' + data_cap_strip[1])) \
                            * (sage.var('s') * sage.var(capacitor_id)) \
                            - sage.var(capacitor_id) \
                            * sage.var('V_initial_' + capacitor_id + '_'
                                 + data_cap_strip[0] + '_'
                                + data_cap_strip[1])
                except KeyError:
                    if data_cap_strip[0] != '0' and data_cap_strip[0] \
                        != 'gnd' and data_cap_strip[0] != 'GND':
                        self.nodal_equations[data_cap_strip[0]] = \
                            (sage.var('V_' + data_cap_strip[0])
                             - sage.var('V_' + data_cap_strip[1])) \
                            * (sage.var('s') * sage.var(capacitor_id)) \
                            - sage.var(capacitor_id) \
                            * sage.var('V_initial_' + capacitor_id + '_'
                                 + data_cap_strip[0] + '_'
                                + data_cap_strip[1])
                try:
                    if data_cap_strip[1] != '0' and data_cap_strip[1] \
                        != 'gnd' and data_cap_strip[1] != 'GND':
                        self.nodal_equations[data_cap_strip[1]] = \
                            self.nodal_equations[data_cap_strip[1]] \
                            + (sage.var('V_' + data_cap_strip[1])
                               - sage.var('V_' + data_cap_strip[0])) \
                            * (sage.var('s') * sage.var(capacitor_id)) \
                            - sage.var(capacitor_id) \
                            * -sage.var('V_initial_' + capacitor_id
                                + '_' + data_cap_strip[0] + '_'
                                + data_cap_strip[1])
                except KeyError:
                    if data_cap_strip[1] != '0' and data_cap_strip[1] \
                        != 'gnd' and data_cap_strip[1] != 'GND':
                        self.nodal_equations[data_cap_strip[1]] = \
                            (sage.var('V_' + data_cap_strip[1])
                             - sage.var('V_' + data_cap_strip[0])) \
                            * (sage.var('s') * sage.var(capacitor_id)) \
                            - sage.var(capacitor_id) \
                            * -sage.var('V_initial_' + capacitor_id
                                + '_' + data_cap_strip[0] + '_'
                                + data_cap_strip[1])

        # set the initial value for the capacitor
        # if set_default_ic_to_zero == True set per default initial value to 0

                if not ignore_all_ic:
                    if set_default_ic_to_zero:

                        self.default_substitutions[sage.var('V_initial_'
                                 + capacitor_id + '_'
                                + data_cap_strip[0] + '_'
                                + data_cap_strip[1])] = '0'
                    try:
                        for key_neq in data_cap_strip[3:]:
                            if key_neq[:3] == 'IC=' or key_neq[:3] \
                                == 'ic=':

                        # set initial condition for the capacitor capacitor_id

                                self.default_substitutions[sage.var('V_initial_'
                                         + capacitor_id + '_'
                                        + data_cap_strip[0] + '_'
                                        + data_cap_strip[1])] = \
                                    key_neq[3:]
                    except:

                # in case of IndexError or other Error do nothing

                        print 'Error in reading initial condition for capacitor ' \
                            + capacitor_id
                        print 'Capacitor data was ',
                        print data_cap_strip
                        raise Exception('Error in reading initial condition for capacitor '
                                 + capacitor_id)
                else:

                    self.default_substitutions[sage.var('V_initial_'
                            + capacitor_id + '_' + data_cap_strip[0]
                            + '_' + data_cap_strip[1])] = '0'
            elif INDUCTOR_EXPR.match(line):

                ind_lineptr = INDUCTOR_EXPR.match(line)
                ind_line_data = line[ind_lineptr.end():]
                inductor_id = line[:ind_lineptr.end()].strip()
                if inductor_id not in self.coupled_inductors_matrix.keys():
                    data_ind = DATA_FIELD.findall(ind_line_data)
                    data_ind_strip = []

                    # data_ind_strip describes the (stripped) data fileds of the 
                    # inductor line line.

                    for key_neq in data_ind:
                        data_ind_strip = data_ind_strip + [key_neq.strip()]

                    # data_ind_strip contains the stripped strings representing 
                    # the data fileds relative to the line line describing 
                    # an inductor in the spice file

                    # the first two data fields describes the nodes to which the 
                    # inductor is connected

                    for key_neq in range(2):
                        if not sage.var('V_' + data_ind_strip[key_neq]) \
                            in self.circuit_variables:
                            self.circuit_variables = self.circuit_variables \
                                + [sage.var('V_' + data_ind_strip[key_neq])]
                            if data_ind_strip[key_neq] != '0' \
                                and data_ind_strip[key_neq] != 'gnd' \
                                and data_ind_strip[key_neq] != 'GND':
                                self.nodal_voltages = self.nodal_voltages \
                                    + [sage.var('V_'
                                       + data_ind_strip[key_neq])]
                    self.circuit_variables = self.circuit_variables \
                        + [sage.var('I_initial_' + inductor_id + '_'
                           + data_ind_strip[0] + '_' + data_ind_strip[1])]
                    self.initial_conditions += [sage.var('I_initial_'
                            + inductor_id + '_' + data_ind_strip[0] + '_'
                            + data_ind_strip[1])]
                    self.circuit_variables = self.circuit_variables \
                        + [sage.var(inductor_id)]
                    try:
                        self.default_substitutions[sage.var(inductor_id)] = \
                            data_ind_strip[2]  
                        self.circuit_graph.add_edge(data_ind_strip[0],
                                data_ind_strip[1], type='L',
                                id=inductor_id, value=data_ind_strip[2:])  
                    except IndexError:
                        self.circuit_graph.add_edge(data_ind_strip[0],
                                data_ind_strip[1], type='L', id=inductor_id)
                    try:
                        if data_ind_strip[0] != '0' and data_ind_strip[0] \
                            != 'gnd' and data_ind_strip[0] != 'GND':
                            self.nodal_equations[data_ind_strip[0]] = \
                                self.nodal_equations[data_ind_strip[0]] \
                                + (sage.var('V_' + data_ind_strip[0])
                                   - sage.var('V_' + data_ind_strip[1])) \
                                / (sage.var('s') * sage.var(inductor_id)) \
                                + sage.var('I_initial_' + inductor_id + '_'
                                    + data_ind_strip[0] + '_'
                                    + data_ind_strip[1]) / sage.var('s')
                    except KeyError:
                        if data_ind_strip[0] != '0' and data_ind_strip[0] \
                            != 'gnd' and data_ind_strip[0] != 'GND':
                            self.nodal_equations[data_ind_strip[0]] = \
                                (sage.var('V_' + data_ind_strip[0])
                                 - sage.var('V_' + data_ind_strip[1])) \
                                / (sage.var('s') * sage.var(inductor_id)) \
                                + sage.var('I_initial_' + inductor_id + '_'
                                    + data_ind_strip[0] + '_'
                                    + data_ind_strip[1]) / sage.var('s')
                    try:
                        if data_ind_strip[1] != '0' and data_ind_strip[1] \
                            != 'gnd' and data_ind_strip[1] != 'GND':
                            self.nodal_equations[data_ind_strip[1]] = \
                                self.nodal_equations[data_ind_strip[1]] \
                                + (sage.var('V_' + data_ind_strip[1])
                                   - sage.var('V_' + data_ind_strip[0])) \
                                / (sage.var('s') * sage.var(inductor_id)) \
                                - sage.var('I_initial_' + inductor_id + '_'
                                    + data_ind_strip[0] + '_'
                                    + data_ind_strip[1]) / sage.var('s')
                    except KeyError:
                        if data_ind_strip[1] != '0' and data_ind_strip[1] \
                            != 'gnd' and data_ind_strip[1] != 'GND':
                            self.nodal_equations[data_ind_strip[1]] = \
                                (sage.var('V_' + data_ind_strip[1])
                                 - sage.var('V_' + data_ind_strip[0])) \
                                / (sage.var('s') * sage.var(inductor_id)) \
                                - sage.var('I_initial_' + inductor_id + '_'
                                    + data_ind_strip[0] + '_'
                                    + data_ind_strip[1]) / sage.var('s')

            # set the initial value for the inductor
            # if set_default_ic_to_zero == True set per default initial value to 0
            # then try: if IC=value parameter is present, adjust initial value 

                    if not ignore_all_ic:
                        if set_default_ic_to_zero:

                            self.default_substitutions[sage.var('I_initial_'
                                     + inductor_id + '_'
                                    + data_ind_strip[0] + '_'
                                    + data_ind_strip[1])] = '0'
                        try:
                            for key_neq in data_ind_strip[3:]:
                                if key_neq[:3] == 'IC=' or key_neq[:3] \
                                    == 'ic=':

                            # set initial condition for the inductor inductor_id

                                    self.default_substitutions[sage.var('I_initial_'
                                             + inductor_id + '_'
                                            + data_ind_strip[0] + '_'
                                            + data_ind_strip[1])] = \
                                        key_neq[3:]
                        except:

                    # in case of IndexError or other Error do nothing

                            print 'Error in reading initial condition for inductor ' \
                                + inductor_id
                            print 'Inductor data was ',
                            print data_ind_strip
                            raise Exception('Error in reading initial condition for inductor '
                                     + inductor_id)
                    else:

                        self.default_substitutions[sage.var('I_initial_'
                                + inductor_id + '_' + data_ind_strip[0]
                                + '_' + data_ind_strip[1])] = '0'

            
            elif VCCS_EXPR.match(line):

        # handle voltage controlled current sources

                vccs_lineptr = VCCS_EXPR.match(line)
                vccs_line_data = line[vccs_lineptr.end():]
                vccs_id = line[:vccs_lineptr.end()].strip()
                data_vccs = DATA_FIELD.findall(vccs_line_data)
                data_vccs_strip = []

        # data_vccs_strip describes the (stripped) data fields of the 
        # voltage controlled current source line line

                for key_neq in data_vccs:
                    data_vccs_strip = data_vccs_strip \
                        + [key_neq.strip()]

        # data_vccs_strip contains the stripped strings representing 
        # the data fields relative to the line line describing a voltage 
        # controlled current source in the spice file
        
        # the first two data fields describes the nodes to which the vccs is 
        # connected. See ngspice documentation

                for key_neq in range(2):
                    if not sage.var('V_' + data_vccs_strip[key_neq]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_'
                               + data_vccs_strip[key_neq])]
                        if data_vccs_strip[key_neq] != '0' \
                            and data_vccs_strip[key_neq] != 'gnd' \
                            and data_vccs_strip[key_neq] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_'
                                   + data_vccs_strip[key_neq])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(vccs_id)]
                try:
                    self.default_substitutions[sage.var(vccs_id)] = \
                        data_vccs_strip[4]
                except IndexError:
                    raise IndexError('Provide a symbolic or numeric value for the vccs '
                             + vccs_id)
                except WrongData, w_dat:
                    raise WrongData('Error in data conversion. Error in converting '
                                     + w_dat.data)

                self.circuit_graph.add_edge(
                    data_vccs_strip[0],
                    data_vccs_strip[1],
                    type='G',
                    id=vccs_id,
                    value=data_vccs_strip[4:],
                    ngspice_line_data=data_vccs_strip[0:4],
                    )
                try:
                    if data_vccs_strip[0] != '0' and data_vccs_strip[0] \
                        != 'gnd' and data_vccs_strip[0] != 'GND':
                        self.nodal_equations[data_vccs_strip[0]] = \
                            self.nodal_equations[data_vccs_strip[0]] \
                            + (sage.var('V_' + data_vccs_strip[2])
                               - sage.var('V_' + data_vccs_strip[3])) \
                            * sage.var(vccs_id)
                except KeyError:

                    if data_vccs_strip[0] != '0' and data_vccs_strip[0] \
                        != 'gnd' and data_vccs_strip[0] != 'GND':
                        self.nodal_equations[data_vccs_strip[0]] = \
                            (sage.var('V_' + data_vccs_strip[2])
                             - sage.var('V_' + data_vccs_strip[3])) \
                            * sage.var(vccs_id)


                try:

                    if data_vccs_strip[1] != '0' and data_vccs_strip[1] \
                        != 'gnd' and data_vccs_strip[1] != 'GND':
                        self.nodal_equations[data_vccs_strip[1]] = \
                            self.nodal_equations[data_vccs_strip[1]] \
                            - (sage.var('V_' + data_vccs_strip[2])
                               - sage.var('V_' + data_vccs_strip[3])) \
                            * sage.var(vccs_id)
                except KeyError:

                    if data_vccs_strip[1] != '0' and data_vccs_strip[1] \
                        != 'gnd' and data_vccs_strip[1] != 'GND':
                        self.nodal_equations[data_vccs_strip[1]] = \
                            -(sage.var('V_' + data_vccs_strip[2])
                              - sage.var('V_' + data_vccs_strip[3])) \
                            * sage.var(vccs_id)

    # following must be at the end of the constructor:
    # finalize equations - must be at the end of the constructor
            

        for key_neq in self.nodal_equations:
            self.nodal_equations[key_neq] = \
                self.nodal_equations[key_neq] == 0

    # substitute the default values in the nodal equations
    # leave this at the end of the constructor -- VERY IMPORTANT !

        self.nodal_equations_substitutions = []
        self.default_substitutions_values = {}
        try:
            for key in self.default_substitutions.keys():
                self.default_substitutions_values[key] = \
                    extract_value(self.default_substitutions[key])

        # set the TEMP temperature and tha values for Kb (Boltzmann constant) 
        # and q (proton charge)
        # do NOT change the order of all these assignments and of this code 
        # (at the end of the constructro).

            self.default_substitutions_values[sage.var('GLOBAL_TEMP'
                    )] = self.temp  # temperature in Kelvin

            self.default_substitutions_values[sage.var('GLOBAL_K_BOLTZMANN'
                    )] = 1.3806488e-23  # Boltzmann constant
            self.default_substitutions_values[sage.var('GLOBAL_q')] = \
                1.602176565e-19  # proton charge....
            self.default_substitutions[sage.var('GLOBAL_TEMP')] = \
                str(self.temp)

            self.default_substitutions[sage.var('GLOBAL_K_BOLTZMANN'
                    )] = '1.3806488e-23'  # Boltzmann constant
            self.default_substitutions[sage.var('GLOBAL_q')] = \
                '1.602176565e-19'  # proton charge....
            self.default_substitutions[sage.var('V_0')] = 0
            self.default_substitutions[sage.var('V_GND')] = 0
            self.default_substitutions[sage.var('V_gnd')] = 0
        except WrongData, w_dat:
            print 'Error in converting value ' \
                + self.default_substitutions[key] + ' associated with ' \
                + key
            print 'Error with ' + w_dat.data
            self.print_information()
            raise WrongData('Error in data conversion. Error with '
                            + w_dat.data)

    # leave this at the very end of the constructor - in the right order !

        subs_zero_ic = {}
        for init_cond in self.initial_conditions:
            try:
                if self.default_substitutions[init_cond] == '0':
                    subs_zero_ic[init_cond] = 0
            except:
                pass

        for key_neq in self.nodal_equations.keys():

            # for loop over the nodal equations

            if not ignore_all_ic:
                self.nodal_equations[key_neq] = \
                    self.nodal_equations[key_neq].substitute({sage.var('V_0'
                        ): 0, sage.var('V_gnd'): 0, sage.var('V_GND'
                        ): 0})
                if set_default_ic_to_zero:
                    self.nodal_equations[key_neq] = \
                        self.nodal_equations[key_neq].substitute(subs_zero_ic)
            else:

            # ignore all Initial Conditions

                subs_zero_ic = {}  # substitutions for zero initial conditions
                for init_cond in self.initial_conditions:
                    subs_zero_ic[init_cond] = 0  # set to zero all initial conditions
                self.nodal_equations[key_neq] = \
                    self.nodal_equations[key_neq].substitute({sage.var('V_0'
                        ): 0, sage.var('V_gnd'): 0, sage.var('V_GND'
                        ): 0}).substitute(subs_zero_ic)  
                        # substitute V_0 and all zero initial conditions 
                        #if ignore_all_ic is true.

        # return numerical nodal equations.

        for key in self.nodal_equations.keys():
            self.nodal_equations_substitutions += \
                [self.nodal_equations[key].substitute(self.additional_equations_explicit).substitute(self.default_substitutions_values)]

        # return a dictonary with the numerical values for 
        # the additional equations

        for key in self.additional_equations_explicit.keys():
            self.additional_equations_exp_sub[key] = \
                self.additional_equations_explicit[key].substitute(self.default_substitutions_values)

        self.circuit_graph_substitutions = self.circuit_graph.copy()
        for key in self.circuit_graph_substitutions.adj.keys():
            for key_2 in \
                self.circuit_graph_substitutions.adj[key].keys():
                for element in \
                    self.circuit_graph_substitutions.adj[key][key_2].keys():
                    if sage.var(self.circuit_graph_substitutions.adj[key][key_2][element]['id'
                                ]) in self.additional_equations_exp_sub:
                        self.circuit_graph_substitutions.adj[key][key_2][element]['value'
                                ] = \
                            [str(self.additional_equations_exp_sub[sage.var(self.circuit_graph_substitutions.adj[key][key_2][element]['id'
                             ])])]

# end of constructor

    def clone(self):
	"""returns a deepcopy of the object
	"""
        return copy.deepcopy(self)

    # when implementing modified nodal analysis, include modified nodal 
    # analysis equations and additional current unknowns

    def solve_nodal_equations_symb(self, set_ind_ss_src_to_zero=False):
        """Returns a symbolic solutions of the nodal equations. 
        
        The solution will be a list containing a list of equations that
        have on the left hand side the nodal voltages and on the right hand 
        side the symbolic expression representing the solution.
        
        set_ind_ss_src_to_zero=True will set all the independent small signal
        sources (in the small signal linearized circuit) to zero when computing
	the solution
        """
        
        if not set_ind_ss_src_to_zero:
            if len(self.nodal_voltages) > 1:
                return sage.solve(self.nodal_equations.values(),
                                  self.nodal_voltages)
            else:
                return [ sage.solve(self.nodal_equations.values(),
                                    self.nodal_voltages) ]
        else:
            nod_eqn = self.nodal_equations.values()
            nod_eqn_adj = []
            for eqn in nod_eqn:
                for i in self.sources_names:
		    #TO DO: CHECK if function or sage.function is correct
                    eqn = eqn.substitute_function(sage.function(i), (0
                            * sage.var('s')).function(sage.var('s')))
                nod_eqn_adj += [eqn]
            if len(self.nodal_voltages) > 1:
                return sage.solve(nod_eqn_adj, self.nodal_voltages)
            else:
                return [ sage.solve(nod_eqn_adj, self.nodal_voltages) ]

    # when implementing modified nodal analysis, include modified nodal analyis
    # equations and additional current unknowns
    

    def solve_nodal_equations_num(self,
            set_ind_ss_src_to_zero=False):
        """Returns a numeric solutions of the nodal equations. 
        
        The solution will be a list containing a list of equations that
        have on the left hand side the nodal voltages and on the right hand 
        side the symbolic expression representing the solution.
        
        set_ind_ss_src_to_zero=True will set all the independent small signal
        sources (in the small signal linearized circuit) to zero when computing
	the solution
        """

        if not set_ind_ss_src_to_zero:
            if len(self.nodal_voltages) > 1:
                return sage.solve(self.nodal_equations_substitutions,
                                  self.nodal_voltages)
            else:
                return [ sage.solve(self.nodal_equations_substitutions,
                                    self.nodal_voltages) ]
        else:
            nod_eqn = self.nodal_equations_substitutions
            nod_eqn_adj = []
            for eqn in nod_eqn:
                for i in self.sources_names:
                    eqn = eqn.substitute_function(sage.function(i), (0
                            * sage.var('s')).function(sage.var('s')))
                nod_eqn_adj += [eqn]
            if len(self.nodal_voltages) > 1:
                return sage.solve(nod_eqn_adj, self.nodal_voltages)
            else:
                return [ sage.solve(nod_eqn_adj, self.nodal_voltages) ]

    # when implementing modified nodal analysis replace all independent voltage
    # sources with a short

        
    def impedance(
        self,
        port,
        with_substitutions=True,
        symbolic=True,
        ):
	"""computes the impedance between two nodes in the circuit

	port=('node1','node2')
	is the pair of nodes considered for the impedance computations, for example
	port=('0','1')
	will compute the impedance between node '1' (if a node '1' is given 
	in the circuit) and ground (node '0').

	with_substitutions=True
	(default) will consider the actual values of all small signal circuit parameters
	that results from the consideration of the operating point of nonlinear devices.
	Otherwise those parameters will be left symbolic.
	(affects the computation of numerical impedance values)

	symbolic=True
	impedance is computed symbolically. (symbolic=False will compute a numeric
	impedance value)

	All internal sources and all the initial conditions will be ignored (set to zero appropriately).

	See also the help for the transimpedance method.
	"""
        try:
            if type(port) != tuple:
                raise Exception("Pass a tuple of two nodes (strings) as an argument in the form ('N1','N2') "
                                )
            elif len(port) != 2:
                raise Exception("Pass a port as a tuple of TWO nodes (string) for calculating an impedance ('N1', 'N2') "
                                )
            elif sage.var('V_' + port[0]) not in self.nodal_voltages \
                + [sage.var('V_0'), sage.var('V_gnd'), sage.var('V_GND'
                   )]:
                raise Exception('First node of the tuple is not a valid node (must be a string describing a node) '
                                )
            elif sage.var('V_' + port[1]) not in self.nodal_voltages \
                + [sage.var('V_0'), sage.var('V_gnd'), sage.var('V_GND'
                   )]:
                raise Exception('Second node of the tuple is not a valid node (must be a string describing a node) '
                                )
        except:

            # re-raise the exception

            raise

        # get the present circuit as a string and remove 
        # the independent sources (for now current sources)
        # when implementing modified nodal analysis replace independent 
        # voltage sources with a short

        circ_without_sources = \
            self.export_lin_circuit_string(remove=self.sources_component_id,
                with_substitutions=with_substitutions)
        if not sage.var('Iss_imp') in self.sources_component_id:
            circ_imp_calc_string = circ_without_sources + 'Iss_imp' \
                + '  ' + port[0] + ' ' + port[1] + ' \n'
        else:
            count = 0
            while sage.var('Iss_imp' + str(count)) \
                in self.sources_component_id:
                count += 1
            circ_imp_calc_string = circ_without_sources + 'Iss_imp' \
                + str(count) + '  ' + port[0] + ' ' + port[1] + ' \n'

        # now we create a circuit for the impedance calculation, with all 
        # initial condition set to zero and containing only the independent 
        # source described with src.

        circ_imp_calc = \
            SmallSignalLinearCircuit(circuit_netlist=circ_imp_calc_string,
                ignore_all_ic=True)
        if symbolic:
            node_voltages = \
                circ_imp_calc.solve_nodal_equations_symb(set_ind_ss_src_to_zero=False)
            try:
                if node_voltages[0][0] == []:
                    raise Exception('An error has occurred')
            except:

                # re-raise exception or index error

                raise
            v_node0 = 0
            v_node1 = 0
            for eqn in node_voltages[0]:
                if sage.var('V_' + port[0]) == eqn.lhs():
                    v_node0 = eqn.rhs()
                if sage.var('V_' + port[1]) == eqn.lhs():
                    v_node1 = eqn.rhs()
            if sage.var('V_' + port[0]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node0 = 0
            if sage.var('V_' + port[1]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node1 = 0
            if len(circ_imp_calc.sources) != 1:
                raise Exception('An error has occurred')
            return (v_node1
                    - v_node0).substitute_function(sage.function(circ_imp_calc.sources_names[0]),
                    (1 + 0 * sage.var('s')).function(sage.var('s')))
        else:
            #calcualte numerical impedance with numerical solutions of the node equations
            node_voltages = \
                circ_imp_calc.solve_nodal_equations_num(set_ind_ss_src_to_zero=False)
            try:
                if node_voltages[0][0] == []:
                    raise Exception('An error has occurred')
            except:

                # re-raise exception or index error

                raise
            v_node0 = 0
            v_node1 = 0
            for eqn in node_voltages[0]:
                if sage.var('V_' + port[0]) == eqn.lhs():
                    v_node0 = eqn.rhs()
                if sage.var('V_' + port[1]) == eqn.lhs():
                    v_node1 = eqn.rhs()
            if sage.var('V_' + port[0]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node0 = 0
            if sage.var('V_' + port[1]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node1 = 0
            if len(circ_imp_calc.sources) != 1:
                raise Exception('An error has occurred')
            return (v_node1
                    - v_node0).substitute_function(sage.function(circ_imp_calc.sources_names[0]),
                    (1 + 0 * sage.var('s')).function(sage.var('s')))
 
            
            
    def transimpedance(self, port_I_in, port_V_out, with_substitutions=True, symbolic=True):
	"""computes the trans-impedance between two nodes in the circuit

	port_I_in=('node1','node2')
	is the pair of nodes considered for the current input
	where the (positive) impressed input current flows into node2 (and out of node1)

	port_V_out=('node3','node4')
	is the pair of nodes considered for the output voltage:
	output_voltage = Vnode4 - Vnode3

	the transimpedance will then be the output voltage/input current

	symbolic=True
	computes a symbolic expression as solution

	with_substitutions=True
	consider the operating point data (affects numerical results with symbolic=False).

	internal sources or initial conditions are ignored.
	"""
    
        try:
            if type(port_I_in) != tuple or type(port_V_out) != tuple:
                raise Exception("Pass a tuple of two nodes (strings) as an argument in the form ('N1','N2') "
                                )
            elif len(port_I_in) != 2 or len(port_V_out) != 2:
                raise Exception("Pass a port as a tuple of TWO nodes (string) for calculating a transimpedance ('N1', 'N2') "
                                )
            elif sage.var('V_' + port_I_in[0]) not in self.nodal_voltages \
                + [sage.var('V_0'), sage.var('V_gnd'), sage.var('V_GND')] or \
                   sage.var('V_' + port_V_out[0]) not in self.nodal_voltages \
                   +[sage.var('V_0'), sage.var('V_gnd'), sage.var('V_GND')]:
                raise Exception('First node of one of the ports (tuples) '
                                'is not a valid node (must be a string describing a node) '
                                )
            elif sage.var('V_' + port_I_in[1]) not in self.nodal_voltages \
                + [sage.var('V_0'), sage.var('V_gnd'), sage.var('V_GND')] or \
                   sage.var('V_' + port_V_out[1]) not in self.nodal_voltages \
                   +[sage.var('V_0'), sage.var('V_gnd'), sage.var('V_GND')]:
                raise Exception('Second node of one of the ports (tuples) '
                                'is not a valid node (must be a string describing a node)'
                                )
        except:

            # re-raise the exception

            raise

        # get the present circuit as a string and remove 
        # the independent sources (for now current sources)
        # when implementing modified nodal analysis replace independent 
        # voltage sources with a short

        circ_without_sources = \
            self.export_lin_circuit_string(remove=self.sources_component_id,
                with_substitutions=with_substitutions)
        if not sage.var('Iss_imp') in self.sources_component_id:
            circ_imp_calc_string = circ_without_sources + 'Iss_imp' \
                + '  ' + port_I_in[0] + ' ' + port_I_in[1] + ' \n'
        else:
            count = 0
            while sage.var('Iss_imp' + str(count)) \
                in self.sources_component_id:
                count += 1
            circ_imp_calc_string = circ_without_sources + 'Iss_imp' \
                + str(count) + '  ' + port_I_in[0] + ' ' + port_I_in[1] + ' \n'

        # now we create a circuit for the transimpedance calculation, with all 
        # initial condition set to zero and containing only the independent 
        # source described with src.

        circ_imp_calc = \
            SmallSignalLinearCircuit(circuit_netlist=circ_imp_calc_string,
                ignore_all_ic=True)
        if symbolic:
            node_voltages = \
                circ_imp_calc.solve_nodal_equations_symb(set_ind_ss_src_to_zero=False)
            try:
                if node_voltages[0][0] == []:
                    raise Exception('An error has occurred')
            except:

                # re-raise exception or index error

                raise
            v_node0 = 0
            v_node1 = 0
            for eqn in node_voltages[0]:
                if sage.var('V_' + port_V_out[0]) == eqn.lhs():
                    v_node0 = eqn.rhs()
                if sage.var('V_' + port_V_out[1]) == eqn.lhs():
                    v_node1 = eqn.rhs()
            if sage.var('V_' + port_V_out[0]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node0 = 0
            if sage.var('V_' + port_V_out[1]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node1 = 0
            if len(circ_imp_calc.sources) != 1:
                raise Exception('An error has occurred')
            return (v_node1
                    - v_node0).substitute_function(sage.function(circ_imp_calc.sources_names[0]),
                    (1 + 0 * sage.var('s')).function(sage.var('s')))
        else:
            #calculate numeric transimpedance with numeric solution of nodal equations 
            node_voltages = \
                circ_imp_calc.solve_nodal_equations_num(set_ind_ss_src_to_zero=False)
            try:
                if node_voltages[0][0] == []:
                    raise Exception('An error has occurred')
            except:

                # re-raise exception or index error

                raise
            v_node0 = 0
            v_node1 = 0
            for eqn in node_voltages[0]:
                if sage.var('V_' + port_V_out[0]) == eqn.lhs():
                    v_node0 = eqn.rhs()
                if sage.var('V_' + port_V_out[1]) == eqn.lhs():
                    v_node1 = eqn.rhs()
            if sage.var('V_' + port_V_out[0]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node0 = 0
            if sage.var('V_' + port_V_out[1]) in [sage.var('V_0'),
                    sage.var('V_gnd'), sage.var('V_GND')]:
                v_node1 = 0
            if len(circ_imp_calc.sources) != 1:
                raise Exception('An error has occurred')
            return (v_node1
                    - v_node0).substitute_function(sage.function(circ_imp_calc.sources_names[0]),
                    (1 + 0 * sage.var('s')).function(sage.var('s')))

    def z_parameters(self, port_in, port_out, with_substitutions=True,
                     symbolic=True):
	"""returns a sage matrix containing the two port network z-parameters
	port_in=('node1', 'node2')
	port_out=('node3', 'node4')

	where nodei is a node identifier

	for any port
	port_x=('nodea', 'nodeb')
	the positive port voltages are nodeb - nodea
	and the positive currents are flowing into nodeb (and out of nodea)
	"""
        z11 = self.impedance(port=port_in, with_substitutions=with_substitutions,
                             symbolic=symbolic)
        z22 = self.impedance(port=port_out, with_substitutions=with_substitutions,
                             symbolic=symbolic)
        z12 = self.transimpedance(port_I_in=port_out, port_V_out=port_in, 
                                  with_substitutions=with_substitutions, 
                                  symbolic=symbolic)
        z21 = self.transimpedance(port_I_in=port_in, port_V_out=port_out, 
                                  with_substitutions=with_substitutions, 
                                  symbolic=symbolic)
        return sage.matrix([[z11, z12],[z21, z22]])


    def dict_default_vals(self):
	"""returns a pair of dictionaries (dict_num, dict_sym)
	where
	dict_num is a dictionary containing as keys the circuit parameters and
	as values the corrisponding numerical values

	dict_sym is a dictionary containing as keys the circuit parameters for
	which a numerical value was not explicitly determined and as values the symbolic
	expressions associated with the considered circuit parameter.

	If all numerical values can be explicitly computed then
	dict_sym will be empty.
	"""
        dict_num = {}
        dict_sym = {}
        for key in self.default_substitutions_values.keys():
            if self.default_substitutions_values[key] in sage.RR:
                dict_num[key] = self.default_substitutions_values[key]
            else:
                if key in self.additional_equations_exp_sub.keys():
                    if self.additional_equations_exp_sub[key] in sage.RR:
                        dict_num[key] = self.additional_equations_exp_sub[key]
                    else:
                        dict_sym[key] = self.additional_equations_exp_sub[key]
                else:
                    dict_sym[key] = self.default_substitutions_values[key]
        return (dict_num, dict_sym)
        
    def print_information(self):
	"""prints some information about the circuit
	"""
        print 'Sources: '
        print self.sources
        print
        print 'Circuit variables: '
        print self.circuit_variables
        print
        print 'Initial conditions: '
        print self.initial_conditions
        print
        print 'Circuit parameters: '
        print self.circuit_parameters
        print
        print 'Default substitutions: '
        print self.default_substitutions
        print
        print 'Nodal voltages: '
        print self.nodal_voltages
        print
        print 'Nodal equations: '
        print self.nodal_equations
        print
        print 'Additional equations: '
        print self.additional_equations
        print
        print 'Additional equations in explicit form: '
        print self.additional_equations_explicit
        print
        print 'Additional equations in explicit form with substitution of default values: '
        print self.additional_equations_exp_sub
        print
        print 'Nodal equations after substitution of default values: '
        print self.nodal_equations_substitutions
        print
        print 'Temperature (Kelvin): '
        print self.temp

    def print_lin_circuit(self):
	"""prints the linear (or linearized) circuit netlist
	"""
        for i in self.circuit_file:
            print i,

    def get_lin_circuit_as_string(self):
	"""returns a string containing the linear circuit 
	(or linearized circuit) netlist
	"""
        ret_val = ''
        for i in self.circuit_file:
            ret_val = ret_val + i
        return ret_val

    def write_lin_circuit_to_file(self, filename):
	"""writes the linear circuit netlist to file 'filename'
	"""
        try:
            lin_cir_file = open(filename, 'w')
        except:
            raise IOError('Cannot write to file ')
        for i in self.circuit_file:
            lin_cir_file.write(i)
        lin_cir_file.close()

    def export_lin_circuit_graph(
        self,
        shorts=None,
        remove=None,
        with_substitutions=True,
        ):
	"""returns a graph (a networkx object) of the linear circuit

	shorts is a dictionary {'node1':'node2', 'node3':'node4'}
	where node1 will be shorted with node2
	and node3 will be shorted with node 4
	etc.
	A key node should not also be a value node for preventing conflicts, so
	when shorting multiple nodes together the proper ordering is important.

	remove is a list of components to remove (replace with open circuit)
	For example remove=['R1','C1']
	will remove the resistor R1 and the capacitor C1
	from the netlist.

	with_substitutions=True
	will explicitly consider numerical values resulting from the operating point
	of nonlinear devices	
	with_substitutions=False will leave the values of parameters resulting from the
	operating point of nonlinear devices as a symbolic expression

	"""
        if shorts == None:
            shorts = {}
        if remove == None:
            remove = []
        try:
            for node in shorts.keys():
                try:
                    edges1 = self.circuit_graph.adj[node][shorts[node]]
                except:
                    edges1 = {}
                try:
                    edges2 = self.circuit_graph.adj[shorts[node]][node]
                except:
                    edges2 = {}
                edges = edges1.values() + edges2.values()
                for edge in edges:
                    try:
                        if edge['type'] == 'L' and edge['coupled_to']['coupled_inductors'].keys() != []:
                            print("Cannot short a coupled inductor")
                            raise WrongUse("WrongUse: Cannot short a coupled inductor")
                    except KeyError:
                        pass
            for edge in remove:
                try:
                    if edge in self.coupled_inductors_matrix.keys():
                        print("Cannot remove a coupled inductor")
                        raise WrongUse("WrongUse: Cannot remove a coupled inductor")
                except IndexError:
                    pass
        except WrongUse:
            raise 
        except:
            raise Exception("An error occurres")
                
        try:
            for node in shorts.keys():
                if not sage.var('V_' + node) in self.nodal_voltages \
                    + [sage.var('V_0'), sage.var('V_gnd'),
                       sage.var('V_GND')]:
                    raise Exception('Dictionary of short circuits not given in a valid form'
                                    )
            for node in shorts.values():
                if not sage.var('V_' + node) in self.nodal_voltages \
                    + [sage.var('V_0'), sage.var('V_gnd'),
                       sage.var('V_GND')]:
                    raise Exception('Dictionary of short circuits not given in a valid form'
                                    )
            for edge in remove:
                if not sage.var(edge) in self.circuit_variables:
                    raise Exception('List of linear components to remove not given in a valid form'
                                    )
            for node in shorts.keys():
                if shorts[node] in shorts.keys():
                    raise Exception('For the dictonary of short circuits a value should not be a key'
                                    )
        except:

            # re-raise the exception

            raise

        try:
            if with_substitutions:
                circ_graph_with_subst = \
                    self.circuit_graph_substitutions.copy()
            else:
                circ_graph_with_subst = self.circuit_graph.copy()
            for edge in circ_graph_with_subst.edges(data=True,
                    keys=True):
                new_node1 = None
                new_node2 = None
                new_attribute_dict = None
                tmp = None
                if edge[0] in shorts.keys():
                    new_node1 = shorts[edge[0]]
                else:
                    new_node1 = edge[0]
                if edge[1] in shorts.keys():
                    new_node2 = shorts[edge[1]]
                else:
                    new_node2 = edge[1]
                new_attribute_dict = edge[3]
                if 'ngspice_line_data' in new_attribute_dict.keys():
                    tmp = []
                    for nod in new_attribute_dict['ngspice_line_data']:
                        if nod in shorts.keys():
                            tmp += [shorts[nod]]
                        else:
                            tmp += [nod]
                    new_attribute_dict['ngspice_line_data'] = tmp
                circ_graph_with_subst.remove_edge(edge[0], edge[1],
                        key=edge[2])
                if not (new_node1 == new_node2
                        and new_attribute_dict['type'] in ['R', 'L', 'C'
                        ]):
                    circ_graph_with_subst.add_edge(new_node1,
                            new_node2, attr_dict=new_attribute_dict)
                else:

                    # else do nothing

                    pass

            for edge in circ_graph_with_subst.edges(data=True,
                    keys=True):
                if edge[3]['id'] in remove:
                    circ_graph_with_subst.remove_edge(edge[0], edge[1],
                            key=edge[2])  # key identifies the edge in a multi-graph
            return circ_graph_with_subst
        except:
            raise Exception('An error occurred when calling export_lin_circuit_graph'
                            )

    def export_lin_circuit_string(
        self,
        shorts=None,
        remove=None,
        with_substitutions=True,
        write_to_file=None,
        ):
	"""returns a string describing the netlist of the linear(ized) circuit

	shorts is a dictionary {'node1':'node2', 'node3':'node4'}
	where node1 will be shorted with node2
	and node3 will be shorted with node 4
	etc.
	A key node should not also be a value node for preventing conflicts, so
	when shorting multiple nodes together the proper ordering is important.

	remove is a list of components to remove (replace with open circuit)
	For example remove=['R1','C1']
	will remove the resistor R1 and the capacitor C1
	from the netlist.

	with_substitutions=True
	will explicitly consider numerical values resulting from the operating point
	of nonlinear devices
	with_substitutions=False will leave the values of parameters resulting from the
	operating point of nonlinear devices as a symbolic expression
	"""

        lin_cir_graph = self.export_lin_circuit_graph(shorts=shorts,
                remove=remove, with_substitutions=with_substitutions)
        circuit_string = ''
        for edge in lin_cir_graph.edges(data=True, keys=True):
            fields_string = ''
            try:
                for field in edge[3]['value']:
                    fields_string += ' ' + field
            except:
                raise Exception('Error in export_lin_circuit_string')
            if not 'ngspice_line_data' in edge[3].keys():
                circuit_string += edge[3]['id'] + ' ' + edge[0] + ' ' \
                    + edge[1] + ' ' + fields_string + ' \n'
            else:
                tmp = ''
                for nod in edge[3]['ngspice_line_data']:
                    tmp += ' ' + nod
                circuit_string += edge[3]['id'] + ' ' + tmp + ' ' \
                    + fields_string + ' \n'
        for key in self.couplings:
            circuit_string += key + ' ' + self.couplings[key]['coupled_inductors'][0].split()[0] \
                    + ' ' + self.couplings[key]['coupled_inductors'][1].split()[0] \
                    + ' ' + self.couplings[key]['coupling_coefficient'] + ' \n'
        if write_to_file != None:
            try:
                out_file = open(write_to_file, 'w')
            except:
                raise IOError('cannot open file')
            out_file.write(circuit_string)
            out_file.close()
        return circuit_string

    def print_symbols(self):
        '''displays information on the symbols adopted by the program (incomplete)'''

        infodict = {
            'V': 'V_nodeid = voltage at node nodeid',
            'V_initial': 'V_initial_capacitorid = initial voltage at capacitor capacitorid according to the ngspice convention',
            'I_initial': 'I_initial_inductorid = initial current in inductor id according to the ngspice convention',
            'BF': 'bjtid_BFop = forward current gain of bjt transistor bjtid',
            }
        keysinfodict = infodict.keys()
        keysinfodict.sort()
        for infokey in keysinfodict:
            print infokey + '  ->  ' + infodict[infokey]

    def write_circuit_graph_dot(
        self,
        output_graph_file='circuit_graph.dot',
        with_substitutions=True,
        edge_labels=True,
        ):
	"""generates a dot file in which a representation of the circuit graph is given.

	output_graoh_file='circuit_graph.dot'
	file in which the graph representation will be written.

	with_substitutions=True
	will consider numerical values resulting from the operating point of nonlinear devices

	edge_labels=True
	will labels edge with informations.
	"""
        try:
            if with_substitutions:
                cir_graph_subst = \
                    sage.Graph(self.circuit_graph_substitutions)
                cir_graph_subst.graphviz_to_file_named(output_graph_file,
                        edge_labels=edge_labels)
            else:
                cir_graph_subst = sage.Graph(self.circuit_graph)
                cir_graph_subst.graphviz_to_file_named(output_graph_file,
                        edge_labels=edge_labels)
        except:
            raise Exception('Error in generating dot file')

    def _preprocessngspicelog(self, logfilelines):
	    position = 0
	    bjt_reached = False
	    for line in logfilelines:
		    position += 1
		    if re.match('^\sBJT: Bipolar Junction Transistor', line):
			    initial_logfile_lines = logfilelines[0:position-1]
			    bjt_reached = True
			    break
	    if bjt_reached:
		    position += 1
		    line = logfilelines[position-1]
		    devices = logfilelines[position-1].split()
		    number_of_devices = len(devices) - 1
		    parameters = []
		    while not re.match('\s*\n', line):
			    position += 1
			    line = logfilelines[position-1]
			    if logfilelines[position-1].split() != []:
				    parameters += [ logfilelines[position-1].split() ]
		    new_data = {} 
		    for bjtnum in xrange(1,number_of_devices+1):
		    	new_data[devices[bjtnum]] = '\tdevice \t' + devices[bjtnum]
		    for parameter in parameters:
			    for bjtnum in xrange(1,number_of_devices+1):
				    new_data[devices[bjtnum]] += '\n'+ \
					'\t' + parameter[0] + '\t' + parameter[bjtnum]
		    return_lines = initial_logfile_lines
		    for bjtnum in xrange(1,number_of_devices+1):
			    return_lines += [' BJT: Bipolar Junction Transistor\n']
			    for newline in new_data[devices[bjtnum]].splitlines():
			    	return_lines += [ newline + '\n' ]
			    return_lines += ['\n']
		    return_lines += logfilelines[position:]
		    return return_lines
	    else:
		    return logfilelines

def extract_value(data):
    '''extracts numerical value from a data field. For example 0.003 from 3m. data must be a stripped string'''

    data_value = \
        re.compile('((\{\w+\})|((\+|-)?\d*\.?\d+(e(\+|-)?\d+))|(\+|-)?\d*\.?\d+(T|G|Meg|K|mil|m|u|n|p|f)?)'
                   )
    data_numfield = re.compile('(\+|-)?\d*\.?\d+')
    if not data_value.match(data):
        raise WrongData(data)
    else:
        if data[0] != '{':
            ptr = data_value.match(data)
            ptrnumfield = data_numfield.match(data)
            stringfactor = data[ptrnumfield.end():ptr.end()]
            if stringfactor != '':
                if stringfactor[0] == 'e':
                    factor = 10 ** int(stringfactor[1:])
                elif stringfactor == 'f':
                    factor = 10 ** -15
                elif stringfactor == 'p':
                    factor = 10 ** -12
                elif stringfactor == 'n':
                    factor = 10 ** -9
                elif stringfactor == 'u':
                    factor = 10 ** -6
                elif stringfactor == 'm':
                    factor = 10 ** -3
                elif stringfactor == 'K':
                    factor = 10 ** 3
                elif stringfactor == 'G':
                    factor = 10 ** 9
                elif stringfactor == 'T':
                    factor = 10 ** 12
                elif stringfactor == 'Meg':
                    factor = 10 ** 6
                elif stringfactor == 'mil':
                    factor = 25.4 * 10 ** -6
                else:
                    factor = 1
            else:
                factor = 1
            return float(data[:ptrnumfield.end()]) * float(factor)
        else:
            return sage.var(data[1:-1])  # intended for symbolic variables

class WrongUse(Exception):
    def __init__(self, data):
        print 'A wrong use was made of the class SmallSignalLinearCircuit'
        print data
        self.data = data
        Exception.__init__(self, data)

class WrongData(Exception):

    def __init__(self, data):
        print 'An exception occurred trying to process the following numerical data value:',
        print data
        self.data = data
        Exception.__init__(self, data)

def simplify_sum(expr, dictionary, treshhold=0):
    """Returns a tuple  (simplified_expr, expr-simplified_expr)
    where simplified_expr is an approximation of expr.	
    
    expr is a symbolic sage expression made of terms added together (a sum of terms)

    dictonary is a dictionary that provides numerical values for the symbolic variables
    in the expression expr

    treshold is a number in [0,1[ that fixes a treshhold.

    If the numerical value of a given term is lower than the numerical value of the term
    that evaluates to the maximum (along all terms) multiplied with the trashhold then the
    corresponding term is ignored.
    """
    if treshhold > 1 or treshhold < 0:
        raise ValueError("treshhold must be between 0 and 1")
    try:
        exp_expr = expr.expand()
        operat = str(exp_expr.operator())
        if operat == '<built-in function add>':
            terms = exp_expr.operands()
            max = 0
            for term in terms:
                if not term.substitute(dictionary) in sage.CC:
                    raise ValueError("Expression expr must be made by terms that when"
                                     " ths substitutions are carried out evaluate"
                                     "to complex numbers. Substitutions are described"
                                     "with the dictionary dictionary.")
                term_value = abs(term.substitute(dictionary))
                if term_value > max:
                    max = term_value
            # max holds the maximim value of the (absolute value of the) 
            # terms in the sum.
            simplified_expr = 0
            for term in terms:
                if abs(term.substitute(dictionary)) < max*treshhold:
                    pass
                else:
                    simplified_expr = simplified_expr + term
        return (simplified_expr, expr-simplified_expr)
                
    except:
        return (expr, 0)
        
def simplify_polynomial(polinom, dictionary, treshhold=0, variable=sage.var('s'), safe_check=True):
    """Returns a tuple  (result, polinom-result)
    where result is an approximation of the polinom polinom given as a sage expression.
    The polynom is considered in the variable variable
    Default: variable=sage.var('s')
    
    dictonary is a dictionary that provides numerical values for the symbolic variables
    in the coefficeints of the polinom.

    treshold is a number in [0,1[ that fixes a treshhold: for treshhold=0 result will be equal
    to polinom.

    Every coefficient of the polynomial (coefficient made by a sum of terms) will be simplified
    invoking simplify_sum.
    """

    if treshhold < 0 or treshhold > 1:
        raise ValueError("treshhold must be between 0 and 1")
    if safe_check:
        if polinom.denominator() != 1:
            raise ValueError("polinom must be a polinom in variable variable")
    try:
        coeff = polinom.coeffs(variable)
        result = 0
        for term in coeff:
            simplified_term = simplify_sum(expr=term[0], dictionary=dictionary, treshhold=treshhold)
            #simplify_sum returns a tuple 
            result = result + (simplified_term[0])*variable**term[1]
        return (result, polinom-result)
    except:
        return (polinom, 0)

def simplify_rational_func(rational_func, dictionary, treshhold=0, variable=sage.var('s'), safe_check=True):
    """simplifies a rational function rational_func (a sage expression)
    in the variable variable (default sage.var('s')).

    dictonary is a dictionary that provides numerical values for the symbolic variables
    in the coefficeints of the polinom.

    treshold is a number in [0,1[ that fixes a treshhold: for treshhold=0 result will be equal
    to polinom.

    A simplified numerator and denominator will be computed invoking simplify_polinomial

    The result will be a dictionaty
    (result, rational_func - result)
    with result a simplified version of rational_func
    """

    if treshhold < 0 or treshhold > 1:
        raise ValueError("treshhold must be between 0 and 1")
    try:
        (num, den) = rational_func.numerator_denominator(safe_check)
        num_simp = simplify_polynomial(polinom=num,dictionary=dictionary,treshhold=treshhold,variable=variable, safe_check=safe_check)
        den_simp = simplify_polynomial(polinom=den,dictionary=dictionary,treshhold=treshhold,variable=variable, safe_check=safe_check)
        # simplify_polynomial retunrs a tuple
        result = num_simp[0]/den_simp[0]
        return (result, rational_func - result)
    except:
        return (rational_func, 0)

def rational_func(large_expression):
    try:
        (num_init_expr, den_init_expr) = large_expression.numerator_denominator(False)
        num_init = num_init_expr.expand()
        den_init = den_init_expr.expand()
        mult = 1
        if str(num_init.operator()) == '<built-in function add>':
            num_init_operands = num_init.operands()
            for term in num_init_operands:
                (term_num, term_denom) = rational_func(term)
                if (mult/term_denom).denominator(True) !=1:
                    mult = mult*term_denom
            if str(den_init.operator()) == '<built-in function add>':
                den_init_operands = den_init.operands()
                for term in den_init_operands:
                    (term_num, term_denom) = rational_func(term)
                    if (mult/term_denom).denominator(True) != 1:
                        mult = mult*term_denom
            else:
                (numerator, denominator) = den_init.numerator_denominator(True)
                if (mult/denominator).denominator(True) != 1:
                    mult = mult*denominator
                den_init_operands = [ den_init ]
            result_numerator = 0
            for term in num_init_operands:
                result_numerator += ((term*mult).expand())
            result_denominator = 0
            for term in den_init_operands:
                result_denominator += ((term*mult).expand())
            return (result_numerator, result_denominator)
        else:
            return large_expression.numerator_denominator(True)
            
    except:
        raise ValueError("Expansion to rational function not successfull")


def get_inductor_data(line, ignore_all_ic, set_default_ic_to_zero):
    return_dict = {}
    return_dict['circuit_variables'] = []
    return_dict['nodal_voltages'] = []
    return_dict['initial_conditions'] = []
    return_dict['default_substitutions'] = {}
    return_dict['line'] = line
    return_dict['value_in_line'] = line.split()[3:]
    lineptr = INDUCTOR_EXPR.match(line)
    line_data = line[lineptr.end():]
    iid = line[:lineptr.end()].strip()
    return_dict['id'] = iid
    return_dict['circuit_variables'] += [ sage.var(iid) ]
    data = DATA_FIELD.findall(line_data)
    data_strip = []
    for element in data:
        data_strip = data_strip + [ element.strip() ]
    return_dict['node0'] = data_strip[0]
    return_dict['node1'] = data_strip[1]
    for index in range(2):
        return_dict['circuit_variables'] += [ sage.var('V_' + data_strip[index]) ]
        if data_strip[index] != '0' and data_strip[index] != 'gnd' and data_strip[index] != 'GND':
            return_dict['nodal_voltages'] += [ sage.var('V_' + data_strip[index]) ]
    return_dict['circuit_variables'] += [ sage.var('I_initial_' + iid + '_' 
                                                   + data_strip[0] + '_' + data_strip[1]) ]
    return_dict['initial_conditions'] += [ sage.var('I_initial_' + iid + '_' 
                                                   + data_strip[0] + '_' + data_strip[1]) ]
    try:
        return_dict['default_substitutions'][sage.var(iid)] = data_strip[2]
        return_dict['value'] = data_strip[2]
    except IndexError:
        pass
    if not ignore_all_ic:
        if set_default_ic_to_zero:
            return_dict['default_substitutions'][sage.var('I_initial_' + iid + '_'
                                                          + data_strip[0] + '_'
                                                          + data_strip[1])] = '0'
        try:
            for field in data_strip[3:]:
                if field[:3] == 'IC=' or field[:3] == 'ic=':
                    return_dict['default_substitutions'][sage.var('I_initial_'
                                                                  + iid + '_'
                                                                  + data_strip[0] + '_'
                                                                  + data_strip[1])] = \
                                                                  field[3:]
            
        except:
            raise Exception('Error in reading initial condition for coupled inductor ' + iid)
    else:    
        return_dict['default_substitutions'][sage.var('I_initial_'
                                                       + iid + '_'
                                                       + data_strip[0] + '_'
                                                       + data_strip[1])] = '0'
                                                       
    return return_dict
