#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Version: 0.2.1.

Author: Alessandro Bernardini

alessandro.bernardini.tum@gmail.com

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

'''

import sage.all as sage
import networkx as nx
import re
import sys
import copy
import string

RESISTOR_EXPR = re.compile('^(\s)*(R|r)\w+')
CAPACITOR_EXPR = re.compile('^(\s)*(C|c)\w+')
INDUCTOR_EXPR = re.compile('^(\s)*(L|l)\w+')
V_EXPR = re.compile('^(\s)*(V|v)\w+')
I_EXPR = re.compile('^(\s)*(I|i)\w+')
VCCS_EXPR = re.compile('^(\s)*(G|g)\w+')
BJT_EXPR = re.compile('^(\s)*(Q|q)\w+')
DATA_FIELD = re.compile('\s+\w+\=-?\w*\.?\w+|\s+-?\w*\.?\w+|\s+\{\w+\}')
TEMP_EXPR = re.compile('^(\s)*(\.temp)\s+\w+', re.IGNORECASE)


class SmallSignalLinearCircuit:

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
        """constructor
    
        with check_operating_region=True a linearized model in dependence of the 
        operating region will be choosen for each semiconductor devices. 
        Otherwise it will be assumed that a standard model is valid. ( for 
        BJT transistors the active region small signal model will be this 
        standard model)
    
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
                    self.spice_batch_output_file = bof.readlines()
                else:
                    self.spice_batch_output_file = []
        else:
            if circuit_netlist == None:
                raise Exception('Provide a circuit file or a linear circuit '
                                'netlist string when constructing a '
                                'SmallSignaLinearCircuit object.'
                                )
            else:

                # circuit_netlist != None

                self.circuit_file = circuit_netlist.splitlines()
                self.original_circuit_file = \
                    circuit_netlist.splitlines()

    # initialize to empty entities and set default temperature.

        self.circuit_graph = nx.MultiDiGraph()
        self.circuit_graph_substitutions = None  # will be generated later
        self.sources = []
        self.sources_names = []
        self.sources_component_id = []
        self.circuit_variables = []
        self.circuit_parameters = []
        self.default_substitutions = {}  # stores the default substitutions with string values
        self.default_substitutions_values = {sage.var('V_0'): 0}  # stores the values for the default substitutions; initialized with V_0 = 0
        self.nodal_equations = {}  # currents flowing out of a node are considered to be positive.
        self.nodal_equations_substitutions = []
        self.additional_equations = []
        self.additional_equations_explicit = {}
        self.additional_equations_exp_sub = {}
        self.initial_conditions = []
        self.thermal_voltages = {}
        self.nodal_voltages = []
        self.operating_regions = {}  # remember the operating regions for the semiconductor devices
        self.temp = 300.15  # default value. Can be overridden if a .temp line is present
        pos = 0

# first for loop over self.circuit_file lines
# in this part preliminary preparatives are taken for a later nodal analysis.

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
                        add_eq_sol = \
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
                        if len(add_eq_sol) == 1:
                            self.additional_equations_explicit.update(add_eq_sol[0])
                        else:

                    # updates the explicit form of the additional equations 
                    # obtained solving for the important quantities

                            raise RuntimeError
                    except:
                        raise Exception('An error has occurred in solving the additional equations for '
                                 + bjt_id)
                else:

                    raise Exception('transistor ' + bjt_id
                                    + ' is not operating in a supported region'
                                    )

# second for  over self.circuit_file lines
# proper nodal analysis done here:

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

                        # set initial condition for the capacitor capacitor_id

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
                    raise Exception('Error in data conversion. Error in converting '
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
            raise Exception('Error in data conversion. Error with '
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
        return copy.deepcopy(self)

    # when implementing modified nodal analysis, include modified nodal 
    # analysis equations and additional current unknowns

    def solve_nodal_equations_symb(self, set_ind_ss_src_to_zero=False):
        """Returns a symbolic solutions of the nodal equations. 
        
        The solution will be a list containing a list of equations that
        have on the left hand side the nodal voltages and on the right hand 
        side the symbolic expression representing the solution.
        
        set_ind_ss_src_to_zero=True will set all the independent small signal
        sources (in the small signal linearized circuit) to zero
        """
        
        if not set_ind_ss_src_to_zero:
            return sage.solve(self.nodal_equations.values(),
                              self.nodal_voltages)
        else:
            nod_eqn = self.nodal_equations.values()
            nod_eqn_adj = []
            for eqn in nod_eqn:
                for i in self.sources_names:
		    #TO DO: CHECK if function or sage.function is correct
                    eqn = eqn.substitute_function(sage.function(i), (0
                            * sage.var('s')).function(sage.var('s')))
                nod_eqn_adj += [eqn]
            return sage.solve(nod_eqn_adj, self.nodal_voltages)

    # when implementing modified nodal analysis, include modified nodal analyis
    # equations and additional current unknowns
    

    def solve_nodal_equations_num(self,
            set_ind_ss_sources_to_zero=False):
        if not set_ind_ss_sources_to_zero:
            return sage.solve(self.nodal_equations_substitutions,
                              self.nodal_voltages)
        else:
            nod_eqn = self.nodal_equations_substitutions
            nod_eqn_adj = []
            for eqn in nod_eqn:
                for i in self.sources_names:
                    eqn = eqn.substitute_function(sage.function(i), (0
                            * sage.var('s')).function(sage.var('s')))
                nod_eqn_adj += [eqn]
            return sage.solve(nod_eqn_adj, self.nodal_voltages)

    # when implementing modified nodal analysis replace all independent voltage
    # sources with a short

        
    def impedance(
        self,
        port,
        with_substitutions=True,
        symbolic=True,
        ):
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
            return (v_node0
                    - v_node1).substitute_function(sage.function(circ_imp_calc.sources_names[0]),
                    (1 + 0 * sage.var('s')).function(sage.var('s')))
        else:

            # to implement

            pass

    def dict_default_vals(self):
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
        for i in self.circuit_file:
            print i,

    def get_lin_circuit_as_string(self):
        ret_val = ''
        for i in self.circuit_file:
            ret_val = ret_val + i
        return ret_val

    def write_lin_circuit_to_file(self, filename):
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
        if shorts == None:
            shorts = {}
        if remove == None:
            remove = []
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
                circuit_string += edge[3]['id'] + '  ' + edge[0] + ' ' \
                    + edge[1] + ' ' + fields_string + ' \n'
            else:
                tmp = ''
                for nod in edge[3]['ngspice_line_data']:
                    tmp += ' ' + nod
                circuit_string += edge[3]['id'] + ' ' + tmp + ' ' \
                    + fields_string + ' \n'
        if write_to_file != None:
            try:
                out_file = open(write_to_file, 'w')
            except:
                raise IOError('cannot open file')
            out_file.write(circuit_string)
            out_file.close()
        return circuit_string

    def print_symbols(self):
        '''displays information on the symbols adopted by the program'''

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


class WrongData:

    def __init__(self, data):
        print 'An exception occurred trying to process the following numerical data value:',
        print data
        self.data = data

def simplify_sum(expr, dict, treshhold=0):
    if treshhold > 1 or treshhold < 0:
        raise ValueError("treshhold must be between 0 and 1")
    try:
        exp_expr = expr.expand()
        operat = str(exp_expr.operator())
        if operat == '<built-in function add>':
            terms = exp_expr.operands()
            max = 0
            for term in terms:
                if not term.substitute(dict) in CC:
                    raise ValueError("Expression expr must be made by terms that when"
                                     " ths substitutions are carried out evaluate"
                                     "to complex numbers. Substitutions are described"
                                     "with the dictionary dict.")
                term_value = abs(term.substitute(dict))
                if term_value > max:
                    max = term_value
            # max holds the maximim value of the (absolute value of the) 
            # terms in the sum.
            simplified_expr = 0
            for term in terms:
                if abs(term.substitute(dict)) < max*treshhold:
                    pass
                else:
                    simplified_expr = simplified_expr + term
        return (simplified_expr, expr-simplified_expr)
                
    except:
        return (expr, 0)
        
def simplify_polynomial(polinom, dict, treshhold=0, variable=sage.var('s'), safe_check=True):
    if treshhold < 0 or treshhold > 1:
        raise ValueError("treshhold must be between 0 and 1")
    if safe_check:
        if polinom.denominator() != 1:
            raise ValueError("polinom must be a polinom in variable variable")
    try:
        coeff = polinom.coeffs(variable)
        result = 0
        for term in coeff:
            simplified_term = simplify_sum(expr=term[0], dict=dict, treshhold=treshhold)
            #simplify_sum returns a tuple 
            result = result + (simplified_term[0])*variable**term[1]
        return (result, polinom-result)
    except:
        return (polinom, 0)

def simplify_rational_func(rational_func, dict, treshhold=0, variable=sage.var('s'), safe_check=True):
    if treshhold < 0 or treshhold > 1:
        raise ValueError("treshhold must be between 0 and 1")
    try:
        (num, den) = rational_func.numerator_denominator()
        num_simp = simplify_polynomial(polinom=num,dict=dict,treshhold=treshhold,variable=variable, safe_check=safe_check)
        den_simp = simplify_polynomial(polinom=den,dict=dict,treshhold=treshhold,variable=variable, safe_check=safe_check)
        # simplify_polynomial retunrs a tuple
        result = num_simp[0]/den_simp[0]
        return (result, rational_func - result)
    except:
        return (rational_func, 0)

        