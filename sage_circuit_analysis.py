#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Version: 0.1 (initial commit).

Author: Alessandro Bernardini

alessandro.bernardini.tum@gmail.com

License: GNU GPL.
http://www.gnu.org/licenses/gpl.html

Disclaimer: THERE IS NO WARRANTY FOR THE PROGRAM (SAPICE and all its provided components), TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS (Alessandro Bernardini) AND/OR OTHER PARTIES PROVIDE THE PROGRAM AS IS WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

THIS PROGRAM WAS NOT WELL TESTED !

Project home:
https://github.com/alessandro-bernardini/SAPICE

Requires: ngspice revision 24; sage version 5.6. Other version of ngspice or sage should work well too.

http://www.sagemath.org/

http://ngspice.sourceforge.net/download.html

See documentation and license

'''

import sage.all as sage
import networkx as nx
import re
import sys

resistor_expr = re.compile('^(\s)*(R|r)\w+')
capacitor_expr = re.compile('^(\s)*(C|c)\w+')
inductor_expr = re.compile('^(\s)*(L|l)\w+')
v_expr = re.compile('^(\s)*(V|v)\w+')
i_expr = re.compile('^(\s)*(I|i)\w+')
vccs_expr = re.compile('^(\s)*(G|g)\w+')
bjt_expr = re.compile('^(\s)*(Q|q)\w+')
data_field = re.compile('\s+\w+\=\w*\.?\w+|\s+\w*\.?\w+|\s+\{\w+\}')
temp_expr = re.compile('^(\s)*(\.temp)\s+\w+', re.IGNORECASE)


class small_signal_linear_circuit:

    def __init__(
        self,
        filename,
        spice_batch_output_file='',
        check_operating_region=True,
        set_default_IC_to_zero=True,
        ignore_all_IC=False,
        ):

    # operating_region is the operating region of the semiconductor devices; with check_operating_region a linearized model in dependence of the operating region will be choosen
    # for each semiconductor devices. Otherwise it will be assumed that a standard model is velid. In this case a BJT transistor will be considered in its active region

        if spice_batch_output_file != '':
            bof = open(spice_batch_output_file, 'r')
            self.spice_batch_output_file = bof.readlines()
        else:
            self.spice_batch_output_file = []
        f = open(filename, 'r')
        self.circuit_file = f.readlines()
	self.original_circuit_file = self.circuit_file
        self.circuit_Graph = nx.MultiGraph()
        self.sources = []
        self.circuit_variables = []
        self.circuit_parameters = []
        self.default_substitutions = {}  # stores the default substitutions with string values
        self.default_substitutions_values = {sage.var('V_0'): '0'}  # stores the values for the default substitutions; initialized with V_0 = 0
        self.nodal_equations = {}  # currents flowing out of a node are considered to be positive.
        self.nodal_equations_substitutions = []
        self.additional_equations = []
        self.additional_equations_explicit = {}
        self.initial_conditions = []
        self.thermal_voltages = {}
        self.nodal_voltages = []
        pos = 0
        self.operating_regions = {}  # remember the operating regions for the semiconductor devices
        self.temp = 300.15  # default value. Can be overridden if a .temp line is present

# first for over self.circuit_file lines
# in this part preliminary preparatives are taken for a later nodal analysis.

        for line in self.circuit_file:
            pos = pos + 1

            if temp_expr.match(line):
                data_temp = data_field.findall(line)  # will include .temp at index 0
                try:
                    temperature = data_temp[0].strip()
                    self.temp = extract_value(temperature) + 273.15  # sets the temperature of the circuit in kelvin
                except:
                    self.temp = 300.15
                    print 'Error in setting the temperature. Default value of 27 C = 300.15 K is assumed'

            if v_expr.match(line):
                print 'Voltage sources are for now not supported. Replace them with an electrically equivalent current source (source transformation)'
                sys.exit()

            if bjt_expr.match(line):
                bjt_lineptr = bjt_expr.match(line)

        # bjt_line_data stores the data of the bjt line

                bjt_line_data = line[bjt_lineptr.end():]
                bjt_id = line[:bjt_lineptr.end()].strip()
                operating_region = 'unknown'

            # if check_operating_region is not set then it will be assumed that the bjt bjt_id operates in its active region

                data_bjt = data_field.findall(bjt_line_data)
                data_bjt_strip = []

                # data_bjt_strip describes the (stripped) data fields of the bjt line line

                for i in data_bjt:
                    data_bjt_strip = data_bjt_strip + [i.strip()]

                # data_bjt_strip contains the stripped strings representing the data fields relative to the line line describing a bjt in the spice file
                # the first 3 data fields describe the nodes to which the bjt is connected

                for i in range(3):
                    if not sage.var('V_' + data_bjt_strip[i]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_bjt_strip[i])]
                        if data_bjt_strip[i] != '0' \
                            and data_bjt_strip[i] != 'gnd' \
                            and data_bjt_strip[i] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_' + data_bjt_strip[i])]

                try:
                    model = data_bjt_strip[3]
                except IndexError:
                    print 'Provide transistor model for transistor connected to nodes ' \
                        + data_bjt_strip[0] + ' ' + data_bjt_strip[1] \
                        + ' ' + data_bjt_strip[2] \
                        + ' (respectively C B E nodes) '
                    sys.exit()
                BJTmodeldescriptionreached = False
                modelforpresentBJTreached = False
                skipnextline = False
                for i in self.spice_batch_output_file:
                    if skipnextline == False:
                        if re.match('^\sBJT\smodels', i) != None:
                            BJTmodeldescriptionreached = True
                            skipnextline = False
                        elif re.match('^\s+model\s+' + model, i,
                                flags=re.IGNORECASE) != None:
                            modelforpresentBJTreached = True
                            skipnextline = True
                        elif BJTmodeldescriptionreached \
                            and modelforpresentBJTreached:
                            matchline = \
                                re.match('^\s+(?P<field>[a-zA-Z]+)\s+(?P<value>(((\+|-)?\d*\.?\d+(e(\+|-)?\d+)?)|(npn)|(pnp)))'
                                    , i, re.IGNORECASE)
                            if matchline != None:
                                if matchline.group('field') != 'type':
                                    self.default_substitutions[bjt_id
        + '_' + matchline.group('field').upper()] = \
    matchline.group('value')
                                else:
                                    if not matchline.group('value') \
    == 'npn':
                                        print 'Only npn bjt are supported '
                                        print matchline.group('field') \
    + ' ' + matchline.group('value') + ' is not supported'
                                        print matchline.group(0)
                                        sys.exit()
                            skipnextline = False
                            if re.match('^\s*\n', i) != None:
                                BJTmodeldescriptionreached = False
                                modelforpresentBJTreached = False
                    else:
                        skipnextline = False
                BJTreached = False
                presentBJTreached = False
                skipnextline = False
                for i in self.spice_batch_output_file:
                    if skipnextline == False:
                        if re.match('^\sBJT:\sBipolar', i) != None:
                            BJTreached = True
                            skipnextline = False
                        elif re.match('^\s+device\s+' + bjt_id, i,
                                flags=re.IGNORECASE) != None:
                            presentBJTreached = True
                            skipnextline = False
                        elif BJTreached and presentBJTreached:
                            matchline = \
                                re.match('^\s+(?P<field>[a-zA-Z]+)\s+(?P<value>(\+|-)?\d*\.?\d+(e(\+|-)?\d+)?)'
                                    , i, re.IGNORECASE)
                            if matchline != None:
                                if matchline.group('field') in ['vbe',
                                        'vbc', 'ic', 'ib']:
                                    self.default_substitutions[bjt_id
        + '_' + matchline.group('field').upper() + '_Q'] = \
    matchline.group('value')
                            skipnextline = False
                        elif re.match('^\s*\n', i) != None:
                            BJTreached = False
                            presentBJTreached = False
                            skipnextline = False
                    else:
                        skipnextline = False
                try:
                    if extract_value(self.default_substitutions[bjt_id
                            + '_VBC_Q']) < 0 \
                        and extract_value(self.default_substitutions[bjt_id
                            + '_VBE_Q']) > 0:
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

            # add equation relating vce with vbe and vbc. This equation is valid in any case, but we still introduce it here in the additional_equations....

                    self.additional_equations = \
                        self.additional_equations + [sage.var(bjt_id
                            + '_VCE_Q') == sage.var(bjt_id + '_VBE_Q')
                            - sage.var(bjt_id + '_VBC_Q')]

            # add equation describing the output resistance of the linearized circuit of the BJT

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

            # add equation describing the capacitance between base and emitter in the linearized circuit

                    self.additional_equations = \
                        self.additional_equations + [sage.var('CBE_'
                            + bjt_id) == sage.var(bjt_id + '_TF')
                            * sage.var('GM_' + bjt_id)]

            # add equation describing the capacitance between base and collector in the linearized circuit

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
                    self.circuit_file = self.circuit_file[:pos-1] + ['RO_'
                             + bjt_id + '  ' + data_bjt_strip[0] + ' '
                            + data_bjt_strip[2] + ' {RO_' + bjt_id
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

                    # updates the explicit form of the additional equations obtained solving for the important quantities

                            raise RuntimeError
                    except:
                        print 'An error has occurred in solving the additional equations for ' \
                            + bjt_id
                        sys.exit()
                else:

                    print 'transistor ' + bjt_id \
                        + ' is not operating in a supported region'
                    sys.exit()

# second for  over self.circuit_file lines
# proper nodal analysis done here:

        for line in self.circuit_file:
            if i_expr.match(line):
                i_lineptr = i_expr.match(line)

        # line describes an independent current source in the spice circuit netlist file

                i_line_data = line[i_lineptr.end():]
                i_id = line[:i_lineptr.end()].strip()
                data_i = data_field.findall(i_line_data)
                data_i_strip = []

        # data_i_strip describes the (stripped) data fields of the current source line line.

                for i in data_i:
                    data_i_strip = data_i_strip + [i.strip()]

        # data_i_strip contains the stripped strings representing the data fields relative to the line line describing a current source in the spice file
        # the first two data fields describes the nodes n+ and n- to which the current source is connected

                for i in range(2):
                    if not sage.var('V_' + data_i_strip[i]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_i_strip[i])]
                        if data_i_strip[i] != '0' and data_i_strip[i] \
                            != 'gnd' and data_i_strip[i] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_' + data_i_strip[i])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(i_id)]

        # only the small signal offset will be considered as a current source. For DC currents (or DC components) a current equal to zero is considered. In case of only a DC current no modification to the nodal equations is needed.

                self.circuit_Graph.add_edge(data_i_strip[0],
                        data_i_strip[1], type='I', id=i_id)
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
        # a current source of positive value forces current to flow out of n+ node, through the current source and into the n- node
        # n+ is associated to dara_i_strip[0] and n- to data_i_strip[1]

                self.sources += [sage.function(i_id + '_small_signal_'
                                 + data_i_strip[0] + '_'
                                 + data_i_strip[1], sage.var('s'))]

        # i_id + '_small_signal_' + data_i_strip[0] + '_' + data_i_strip[1] is the current imposed by source i_id flowing from data_i_strip[0] to data_i_strip[1] through the current source i_id itself.........

            if resistor_expr.match(line):
                res_lineptr = resistor_expr.match(line)

                # line describes a resistor in the spice circuit netlist file

                res_line_data = line[res_lineptr.end():]
                resistor_id = line[:res_lineptr.end()].strip()
                data_res = data_field.findall(res_line_data)
                data_res_strip = []

                # data_res_strip describes the (stripped) data fileds of the resistor line line.

                for i in data_res:
                    data_res_strip = data_res_strip + [i.strip()]

                # data_res_strip contains the stripped strings representing the data fileds relative to the line line describing a resistor in the spice file
                # the first two data fields describes the nodes to which the resistor is connected

                for i in range(2):
                    if not sage.var('V_' + data_res_strip[i]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_res_strip[i])]
                        if data_res_strip[i] != '0' \
                            and data_res_strip[i] != 'gnd' \
                            and data_res_strip[i] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_' + data_res_strip[i])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(resistor_id)]
                try:
                    self.default_substitutions[resistor_id] = \
                        data_res_strip[2]  # adjust this for extracting numerical value
                    self.circuit_Graph.add_edge(data_res_strip[0],
                            data_res_strip[1], type='R',
                            id=resistor_id, value=data_res_strip[2:])  # adjust this for extracting numerical value
                except IndexError:
                    self.circuit_Graph.add_edge(data_res_strip[0],
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
            elif capacitor_expr.match(line):
                cap_lineptr = capacitor_expr.match(line)
                cap_line_data = line[cap_lineptr.end():]
                capacitor_id = line[:cap_lineptr.end()].strip()
                data_cap = data_field.findall(cap_line_data)
                data_cap_strip = []

                # data_cap_strip describes the (stripped) data fileds of the capacitor line line.

                for i in data_cap:
                    data_cap_strip = data_cap_strip + [i.strip()]

                # data_cap_strip contains the stripped strings representing the data fileds relative to the line line describing a capacitor in the spice file
                # the first two data fields describes the nodes to which the capacitor is connected

                for i in range(2):
                    if not sage.var('V_' + data_cap_strip[i]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_cap_strip[i])]
                        if data_cap_strip[i] != '0' \
                            and data_cap_strip[i] != 'gnd' \
                            and data_cap_strip[i] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_' + data_cap_strip[i])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var('V_initial_' + capacitor_id + '_'
                       + data_cap_strip[0] + '_' + data_cap_strip[1])]
                self.initial_conditions += [sage.var('V_initial_'
                        + capacitor_id + '_' + data_cap_strip[0] + '_'
                        + data_cap_strip[1])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(capacitor_id)]
                try:
                    self.default_substitutions[capacitor_id] = \
                        data_cap_strip[2]  # adjust this for extracting numerical value
                    self.circuit_Graph.add_edge(data_cap_strip[0],
                            data_cap_strip[1], type='C',
                            id=capacitor_id, value=data_cap_strip[2:])  # adjust this for extracting numerical value
                except IndexError:
                    self.circuit_Graph.add_edge(data_cap_strip[0],
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
                            * (sage.var('s') * sage.var(resistor_id)) \
                            - sage.var(capacitor_id) \
                            * -sage.var('V_initial_' + capacitor_id
                                + '_' + data_cap_strip[0] + '_'
                                + data_cap_strip[1])

        # set the initial value for the capacitor
        # if set_default_IC_to_zero == True set per default initial value to zero
        # then try: if IC=value parameter is present, adjust initial value to given value

                if not ignore_all_IC:
                    if set_default_IC_to_zero:
                        self.default_substitutions[sage.var('V_initial_'
                                 + capacitor_id + '_'
                                + data_cap_strip[0] + '_'
                                + data_cap_strip[1])] = '0'
                    try:
                        for i in data_cap_strip[3:]:
                            if i[:3] == 'IC=' or i[:3] == 'ic=':

                        # set initial condition for the capacitor capacitor_id

                                self.default_substitutions[sage.var('V_initial_'
                                         + capacitor_id + '_'
                                        + data_cap_strip[0] + '_'
                                        + data_cap_strip[1])] = i[3:]
                    except:

                # in case of IndexError or other Error do nothing

                        print 'Error in reading initial condition for capacitor ' \
                            + capacitor_id
                        print 'Capacitor data was ',
                        print data_cap_strip
                        sys.exit()
                else:
                    self.default_substitutions[sage.var('V_initial_'
                            + capacitor_id + '_' + data_cap_strip[0]
                            + '_' + data_cap_strip[1])] = '0'
            elif inductor_expr.match(line):

                ind_lineptr = inductor_expr.match(line)
                ind_line_data = line[ind_lineptr.end():]
                inductor_id = line[:ind_lineptr.end()].strip()
                data_ind = data_field.findall(ind_line_data)
                data_ind_strip = []

                # data_ind_strip describes the (stripped) data fileds of the inductor line line.

                for i in data_ind:
                    data_ind_strip = data_ind_strip + [i.strip()]

                # data_ind_strip contains the stripped strings representing the data fileds relative to the line line describing an inductor in the spice file
                # the first two data fields describes the nodes to which the inductor is connected

                for i in range(2):
                    if not sage.var('V_' + data_ind_strip[i]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_ind_strip[i])]
                        if data_ind_strip[i] != '0' \
                            and data_ind_strip[i] != 'gnd' \
                            and data_ind_strip[i] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_' + data_ind_strip[i])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var('I_initial_' + inductor_id + '_'
                       + data_ind_strip[0] + '_' + data_ind_strip[1])]
                self.initial_conditions += [sage.var('I_initial_'
                        + inductor_id + '_' + data_ind_strip[0] + '_'
                        + data_ind_strip[1])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(inductor_id)]
                try:
                    self.default_substitutions[inductor_id] = \
                        data_ind_strip[2]  # adjust this for extracting numerical value
                    self.circuit_Graph.add_edge(data_ind_strip[0],
                            data_ind_strip[1], type='L',
                            id=inductor_id, value=data_ind_strip[2:])  # adjust this for extracting numerical value
                except IndexError:
                    self.circuit_Graph.add_edge(data_ind_strip[0],
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
        # if set_default_IC_to_zero == True set per default initial value to zero
        # then try: if IC=value parameter is present, adjust initial value to given value

                if not ignore_all_IC:
                    if set_default_IC_to_zero:
                        self.default_substitutions[sage.var('I_initial_'
                                 + inductor_id + '_'
                                + data_ind_strip[0] + '_'
                                + data_ind_strip[1])] = '0'
                    try:
                        for i in data_ind_strip[3:]:
                            if i[:3] == 'IC=' or i[:3] == 'ic=':

                        # set initial condition for the capacitor capacitor_id

                                self.default_substitutions[sage.var('I_initial_'
                                         + inductor_id + '_'
                                        + data_ind_strip[0] + '_'
                                        + data_ind_strip[1])] = i[3:]
                    except:

                # in case of IndexError or other Error do nothing

                        print 'Error in reading initial condition for inductor ' \
                            + inductor_id
                        print 'Inductor data was ',
                        print data_ind_strip
                        sys.exit()
                else:
                    self.default_substitutions[sage.var('I_initial_'
                            + inductor_id + '_' + data_ind_strip[0]
                            + '_' + data_ind_strip[1])] = '0'
            elif vccs_expr.match(line):

        # handle voltage controlled current sources

                vccs_lineptr = vccs_expr.match(line)
                vccs_line_data = line[vccs_lineptr.end():]
                vccs_id = line[:vccs_lineptr.end()].strip()
                data_vccs = data_field.findall(vccs_line_data)
                data_vccs_strip = []

        # data_vccs_strip describes the (stripped) data fields of the voltage controlled current source line line

                for i in data_vccs:
                    data_vccs_strip = data_vccs_strip + [i.strip()]

        # data_vccs_strip contains the stripped strings representing the data fields relative to the line line describing a voltage controlled current source in the spice file
        # the first two data fields describes the nodes to which the vccs is connected. See ngspice documentation

                for i in range(2):
                    if not sage.var('V_' + data_vccs_strip[i]) \
                        in self.circuit_variables:
                        self.circuit_variables = self.circuit_variables \
                            + [sage.var('V_' + data_vccs_strip[i])]
                        if data_vccs_strip[i] != '0' \
                            and data_vccs_strip[i] != 'gnd' \
                            and data_vccs_strip[i] != 'GND':
                            self.nodal_voltages = self.nodal_voltages \
                                + [sage.var('V_' + data_vccs_strip[i])]
                self.circuit_variables = self.circuit_variables \
                    + [sage.var(vccs_id)]
                self.circuit_Graph.add_edge(data_vccs_strip[1],
                        data_vccs_strip[0], type='G', id=vccs_id,
                        value=data_vccs_strip[2:])  # adjust this for extracting numerical value
                try:
                    if data_vccs_strip[0] != '0' and data_vccs_strip[0] \
                        != 'gnd' and data_vccs_strip[0] != 'GND':
                        self.nodal_equations[data_vccs_strip[0]] = \
                            self.nodal_equations[data_vccs_strip[0]] \
                            + (sage.var('V_' + data_vccs_strip[2])
                               - sage.var('V_' + data_vccs_strip[3])) \
                            * extract_value(data_vccs_strip[4])
                except KeyError:

                    if data_vccs_strip[0] != '0' and data_vccs_strip[0] \
                        != 'gnd' and data_vccs_strip[0] != 'GND':
                        self.nodal_equations[data_vccs_strip[0]] = \
                            (sage.var('V_' + data_vccs_strip[2])
                             - sage.var('V_' + data_vccs_strip[3])) \
                            * extract_value(data_vccs_strip[4])
                except WrongData:

                    print 'Error with the value of transconductance of voltage controlled current source ' \
                        + vccs_id
                    print 'Data was '
                    print data_vccs_strip
                    sys.exit()

                try:

                    if data_vccs_strip[1] != '0' and data_vccs_strip[1] \
                        != 'gnd' and data_vccs_strip[1] != 'GND':
                        self.nodal_equations[data_vccs_strip[1]] = \
                            self.nodal_equations[data_vccs_strip[1]] \
                            - (sage.var('V_' + data_vccs_strip[2])
                               - sage.var('V_' + data_vccs_strip[3])) \
                            * extract_value(data_vccs_strip[4])
                except KeyError:
                    if data_vccs_strip[1] != '0' and data_vccs_strip[1] \
                        != 'gnd' and data_vccs_strip[1] != 'GND':
                        self.nodal_equations[data_vccs_strip[1]] = \
                            -(sage.var('V_' + data_vccs_strip[2])
                              - sage.var('V_' + data_vccs_strip[3])) \
                            * extract_value(data_vccs_strip[4])
                except WrongData:
                    print 'Error with the value of transconductance of voltage controlled current source ' \
                        + vccs_id
                    print 'Data was '
                    print data_vccs_strip
                    sys.exit()

    # following must be at the end of the constructor:
    # finalize equations - must be at the end of the constructor

        for i in self.nodal_equations:
            self.nodal_equations[i] = self.nodal_equations[i] == 0

    # substitute the default values in the nodal equations
    # leave this at the end of the constructor -- VERY IMPORTANT !

        self.nodal_equations_substitutions = []
        self.default_substitutions_values = {}
        try:
            for key in self.default_substitutions.keys():
                self.default_substitutions_values[sage.var(key)] = \
                    extract_value(self.default_substitutions[key])

        # set the TEMP temperature and tha values for Kb (Boltzmann constant) and q (proton charge)
        # do NOT change the order of all these assignments and of this code (at the end of the constructro).

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
        except WrongData, wd:
            print 'Error in converting value ' \
                + self.default_substitutions[key] + ' associated with ' \
                + key
            print 'Error with ' + wd.data
            self.print_information()
            sys.exit()

    # leave this at the very end of the constructor - in the right order !

        for i in self.nodal_equations.keys():
            if not ignore_all_IC:
                self.nodal_equations[i] = \
                    self.nodal_equations[i].substitute({sage.var('V_0'
                        ): 0, sage.var('V_gnd'): 0, sage.var('V_GND'
                        ): 0})
            else:

            # ignore all Initial Conditions

                subs_zero_IC = {}  # substitutions for zero initial conditions
                for i in self.initial_conditions:
                    subs_zero_IC[i] = 0
                for key in self.nodal_equations.keys():
                    self.nodal_equations[key] = \
                        self.nodal_equations[key].substitute({sage.var('V_0'
                            ): 0, sage.var('V_gnd'): 0, sage.var('V_GND'
                            ): 0}).substitute(subs_zero_IC)  # substitute V_0 and the zero initial conditions if ignora_all_IC is true.

        for key in self.nodal_equations.keys():
            self.nodal_equations_substitutions += \
                [self.nodal_equations[key].substitute(self.additional_equations_explicit).substitute(self.default_substitutions_values)]

        # first substitute additional_equations_explicit and then the numerical values

    def solve_nodal_equations(self):
        return sage.solve(self.nodal_equations.values(),
                          self.nodal_voltages)

    def solve_nodal_equations_with_substitutions(self):
        return sage.solve(self.nodal_equations_substitutions,
                          self.nodal_voltages)

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
        print 'Nodal equations after substitution of default values: '
        print self.nodal_equations_substitutions
        print
        print 'Temperature (Kelvin): '
        print self.temp

    def print_lin_circuit(self):
        for i in self.circuit_file:
            print i,

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

    def create_circuit_graph(self, output_graph_file='circuit_graph.dot'
                             ):
        pass


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


# to do: implement voltage sources and current sources

