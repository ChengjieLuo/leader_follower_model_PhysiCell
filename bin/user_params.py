 
# This file is auto-generated from a Python script that parses a PhysiCell configuration (.xml) file.
#
# Edit at your own risk.
#
import os
from ipywidgets import Label,Text,Checkbox,Button,HBox,VBox,FloatText,IntText,BoundedIntText,BoundedFloatText,Layout,Box
    
class UserTab(object):

    def __init__(self):
        
        micron_units = Label('micron')   # use "option m" (Mac, for micro symbol)

        constWidth = '180px'
        tab_height = '500px'
        stepsize = 10

        #style = {'description_width': '250px'}
        style = {'description_width': '25%'}
        layout = {'width': '400px'}

        name_button_layout={'width':'25%'}
        widget_layout = {'width': '15%'}
        units_button_layout ={'width':'15%'}
        desc_button_layout={'width':'45%'}

        param_name1 = Button(description='Model_Initialization', disabled=True, layout=name_button_layout)
        param_name1.style.button_color = 'lightgreen'

        self.Model_Initialization = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name2 = Button(description='ECM_initialization_parameters', disabled=True, layout=name_button_layout)
        param_name2.style.button_color = 'tan'

        self.ECM_initialization_parameters = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name3 = Button(description='initial_anisotropy', disabled=True, layout=name_button_layout)
        param_name3.style.button_color = 'lightgreen'

        self.initial_anisotropy = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name4 = Button(description='initial_ECM_density', disabled=True, layout=name_button_layout)
        param_name4.style.button_color = 'tan'

        self.initial_ECM_density = FloatText(
          value=0.5,
          step=0.1,
          style=style, layout=widget_layout)

        param_name5 = Button(description='Lesion_intialization', disabled=True, layout=name_button_layout)
        param_name5.style.button_color = 'lightgreen'

        self.Lesion_intialization = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name6 = Button(description='tumor_radius', disabled=True, layout=name_button_layout)
        param_name6.style.button_color = 'tan'

        self.tumor_radius = FloatText(
          value=175,
          step=10,
          style=style, layout=widget_layout)

        param_name7 = Button(description='initial_leader_cell_fraction', disabled=True, layout=name_button_layout)
        param_name7.style.button_color = 'lightgreen'

        self.initial_leader_cell_fraction = FloatText(
          value=0.2,
          step=0.01,
          style=style, layout=widget_layout)

        param_name8 = Button(description='Cell_properties', disabled=True, layout=name_button_layout)
        param_name8.style.button_color = 'tan'

        self.Cell_properties = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name9 = Button(description='oxygen_uptake', disabled=True, layout=name_button_layout)
        param_name9.style.button_color = 'lightgreen'

        self.oxygen_uptake = FloatText(
          value=10,
          step=1,
          style=style, layout=widget_layout)

        param_name10 = Button(description='leader_adhesion', disabled=True, layout=name_button_layout)
        param_name10.style.button_color = 'tan'

        self.leader_adhesion = FloatText(
          value=5.0,
          step=0.1,
          style=style, layout=widget_layout)

        param_name11 = Button(description='leader_repulsion', disabled=True, layout=name_button_layout)
        param_name11.style.button_color = 'lightgreen'

        self.leader_repulsion = FloatText(
          value=25,
          step=1,
          style=style, layout=widget_layout)

        param_name12 = Button(description='follower_adhesion', disabled=True, layout=name_button_layout)
        param_name12.style.button_color = 'tan'

        self.follower_adhesion = FloatText(
          value=5.0,
          step=0.1,
          style=style, layout=widget_layout)

        param_name13 = Button(description='follower_repulsion', disabled=True, layout=name_button_layout)
        param_name13.style.button_color = 'lightgreen'

        self.follower_repulsion = FloatText(
          value=25,
          step=1,
          style=style, layout=widget_layout)

        param_name14 = Button(description='default_cell_speed', disabled=True, layout=name_button_layout)
        param_name14.style.button_color = 'tan'

        self.default_cell_speed = FloatText(
          value=0.5,
          step=0.1,
          style=style, layout=widget_layout)

        param_name15 = Button(description='rho_L', disabled=True, layout=name_button_layout)
        param_name15.style.button_color = 'lightgreen'

        self.rho_L = FloatText(
          value=0.0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name16 = Button(description='rho_H', disabled=True, layout=name_button_layout)
        param_name16.style.button_color = 'tan'

        self.rho_H = FloatText(
          value=1.0,
          step=0.1,
          style=style, layout=widget_layout)

        param_name17 = Button(description='rho_I', disabled=True, layout=name_button_layout)
        param_name17.style.button_color = 'lightgreen'

        self.rho_I = FloatText(
          value=0.5,
          step=0.1,
          style=style, layout=widget_layout)

        param_name18 = Button(description='default_persistence_time', disabled=True, layout=name_button_layout)
        param_name18.style.button_color = 'tan'

        self.default_persistence_time = FloatText(
          value=10.0,
          step=1,
          style=style, layout=widget_layout)

        param_name19 = Button(description='default_ECM_density_target', disabled=True, layout=name_button_layout)
        param_name19.style.button_color = 'lightgreen'

        self.default_ECM_density_target = FloatText(
          value=0.5,
          step=0.1,
          style=style, layout=widget_layout)

        param_name20 = Button(description='default_ECM_sensitivity', disabled=True, layout=name_button_layout)
        param_name20.style.button_color = 'tan'

        self.default_ECM_sensitivity = FloatText(
          value=1.0,
          step=0.1,
          style=style, layout=widget_layout)

        param_name21 = Button(description='default_hysteresis_bias', disabled=True, layout=name_button_layout)
        param_name21.style.button_color = 'lightgreen'

        self.default_hysteresis_bias = FloatText(
          value=1.0,
          step=0.1,
          style=style, layout=widget_layout)

        param_name22 = Button(description='oxygen_migration_bias_for_leaders', disabled=True, layout=name_button_layout)
        param_name22.style.button_color = 'tan'

        self.oxygen_migration_bias_for_leaders = FloatText(
          value=0.95,
          step=0.1,
          style=style, layout=widget_layout)

        param_name23 = Button(description='Microenvironment_parameters', disabled=True, layout=name_button_layout)
        param_name23.style.button_color = 'lightgreen'

        self.Microenvironment_parameters = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name24 = Button(description='chemotactic_substrate_decay_rate', disabled=True, layout=name_button_layout)
        param_name24.style.button_color = 'tan'

        self.chemotactic_substrate_decay_rate = FloatText(
          value=0.1,
          step=0.01,
          style=style, layout=widget_layout)

        param_name25 = Button(description='Model_Selection_Parameters', disabled=True, layout=name_button_layout)
        param_name25.style.button_color = 'lightgreen'

        self.Model_Selection_Parameters = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name26 = Button(description='cell_motility_ECM_interaction_model_selector', disabled=True, layout=name_button_layout)
        param_name26.style.button_color = 'tan'

        self.cell_motility_ECM_interaction_model_selector = Text(
          value='previous motility vector based memory',
          style=style, layout=widget_layout)

        param_name27 = Button(description='Testing_and_setup_parameters', disabled=True, layout=name_button_layout)
        param_name27.style.button_color = 'lightgreen'

        self.Testing_and_setup_parameters = FloatText(
          value=0,
          step=0.01,
          style=style, layout=widget_layout)

        param_name28 = Button(description='normalize_ECM_influenced_motility_vector', disabled=True, layout=name_button_layout)
        param_name28.style.button_color = 'tan'

        self.normalize_ECM_influenced_motility_vector = Checkbox(
          value=False,
          style=style, layout=widget_layout)

        param_name29 = Button(description='duration_of_uE_conditioning', disabled=True, layout=name_button_layout)
        param_name29.style.button_color = 'lightgreen'

        self.duration_of_uE_conditioning = FloatText(
          value=10,
          step=1,
          style=style, layout=widget_layout)

        param_name30 = Button(description='freeze_uE_profile', disabled=True, layout=name_button_layout)
        param_name30.style.button_color = 'tan'

        self.freeze_uE_profile = Checkbox(
          value=False,
          style=style, layout=widget_layout)

        param_name31 = Button(description='unit_test_setup', disabled=True, layout=name_button_layout)
        param_name31.style.button_color = 'lightgreen'

        self.unit_test_setup = IntText(
          value=0,
          step=1,
          style=style, layout=widget_layout)

        param_name32 = Button(description='march_unit_test_setup', disabled=True, layout=name_button_layout)
        param_name32.style.button_color = 'tan'

        self.march_unit_test_setup = IntText(
          value=0,
          step=1,
          style=style, layout=widget_layout)

        param_name33 = Button(description='cell_setup', disabled=True, layout=name_button_layout)
        param_name33.style.button_color = 'lightgreen'

        self.cell_setup = Text(
          value='lesion',
          style=style, layout=widget_layout)

        param_name34 = Button(description='ECM_orientation_setup', disabled=True, layout=name_button_layout)
        param_name34.style.button_color = 'tan'

        self.ECM_orientation_setup = Text(
          value='random',
          style=style, layout=widget_layout)

        param_name35 = Button(description='chemical_field_setup', disabled=True, layout=name_button_layout)
        param_name35.style.button_color = 'lightgreen'

        self.chemical_field_setup = Text(
          value='none',
          style=style, layout=widget_layout)

        param_name36 = Button(description='angle_of_chemical_field_gradient', disabled=True, layout=name_button_layout)
        param_name36.style.button_color = 'tan'

        self.angle_of_chemical_field_gradient = FloatText(
          value=45.0,
          step=1,
          style=style, layout=widget_layout)

        units_button1 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button1.style.button_color = 'lightgreen'
        units_button2 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button2.style.button_color = 'tan'
        units_button3 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button3.style.button_color = 'lightgreen'
        units_button4 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button4.style.button_color = 'tan'
        units_button5 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button5.style.button_color = 'lightgreen'
        units_button6 = Button(description='microns', disabled=True, layout=units_button_layout) 
        units_button6.style.button_color = 'tan'
        units_button7 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button7.style.button_color = 'lightgreen'
        units_button8 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button8.style.button_color = 'tan'
        units_button9 = Button(description='mmHg/min', disabled=True, layout=units_button_layout) 
        units_button9.style.button_color = 'lightgreen'
        units_button10 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button10.style.button_color = 'tan'
        units_button11 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button11.style.button_color = 'lightgreen'
        units_button12 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button12.style.button_color = 'tan'
        units_button13 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button13.style.button_color = 'lightgreen'
        units_button14 = Button(description='microns/minute', disabled=True, layout=units_button_layout) 
        units_button14.style.button_color = 'tan'
        units_button15 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button15.style.button_color = 'lightgreen'
        units_button16 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button16.style.button_color = 'tan'
        units_button17 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button17.style.button_color = 'lightgreen'
        units_button18 = Button(description='minutes', disabled=True, layout=units_button_layout) 
        units_button18.style.button_color = 'tan'
        units_button19 = Button(description='dimenionless', disabled=True, layout=units_button_layout) 
        units_button19.style.button_color = 'lightgreen'
        units_button20 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button20.style.button_color = 'tan'
        units_button21 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button21.style.button_color = 'lightgreen'
        units_button22 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button22.style.button_color = 'tan'
        units_button23 = Button(description='dimensionlesss', disabled=True, layout=units_button_layout) 
        units_button23.style.button_color = 'lightgreen'
        units_button24 = Button(description='1/minutes', disabled=True, layout=units_button_layout) 
        units_button24.style.button_color = 'tan'
        units_button25 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button25.style.button_color = 'lightgreen'
        units_button26 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button26.style.button_color = 'tan'
        units_button27 = Button(description='dimensionlesss', disabled=True, layout=units_button_layout) 
        units_button27.style.button_color = 'lightgreen'
        units_button28 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button28.style.button_color = 'tan'
        units_button29 = Button(description='minutes', disabled=True, layout=units_button_layout) 
        units_button29.style.button_color = 'lightgreen'
        units_button30 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button30.style.button_color = 'tan'
        units_button31 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button31.style.button_color = 'lightgreen'
        units_button32 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button32.style.button_color = 'tan'
        units_button33 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button33.style.button_color = 'lightgreen'
        units_button34 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button34.style.button_color = 'tan'
        units_button35 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button35.style.button_color = 'lightgreen'
        units_button36 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button36.style.button_color = 'tan'

        desc_button1 = Button(description='This is a spacer', disabled=True, layout=desc_button_layout) 
        desc_button1.style.button_color = 'lightgreen'
        desc_button2 = Button(description='ECM_initialization_parameters', disabled=True, layout=desc_button_layout) 
        desc_button2.style.button_color = 'tan'
        desc_button3 = Button(description='Initial ECM anisotropy', disabled=True, layout=desc_button_layout) 
        desc_button3.style.button_color = 'lightgreen'
        desc_button4 = Button(description='Initial ECM density', disabled=True, layout=desc_button_layout) 
        desc_button4.style.button_color = 'tan'
        desc_button5 = Button(description='Lesion_intialization', disabled=True, layout=desc_button_layout) 
        desc_button5.style.button_color = 'lightgreen'
        desc_button6 = Button(description='Initial cell mass radius', disabled=True, layout=desc_button_layout) 
        desc_button6.style.button_color = 'tan'
        desc_button7 = Button(description='Average percentage of leaders initialized', disabled=True, layout=desc_button_layout) 
        desc_button7.style.button_color = 'lightgreen'
        desc_button8 = Button(description='Cell_properties', disabled=True, layout=desc_button_layout) 
        desc_button8.style.button_color = 'tan'
        desc_button9 = Button(description='Default cell oxygen uptake rate', disabled=True, layout=desc_button_layout) 
        desc_button9.style.button_color = 'lightgreen'
        desc_button10 = Button(description='Strength of leader cells\' cell-cell adhesion', disabled=True, layout=desc_button_layout) 
        desc_button10.style.button_color = 'tan'
        desc_button11 = Button(description='Strength of leader cells\' cell-cell repulsion', disabled=True, layout=desc_button_layout) 
        desc_button11.style.button_color = 'lightgreen'
        desc_button12 = Button(description='Strength of follower cells\' cell-cell adhesion', disabled=True, layout=desc_button_layout) 
        desc_button12.style.button_color = 'tan'
        desc_button13 = Button(description='Strength of follower cells\' cell-cell repulsion', disabled=True, layout=desc_button_layout) 
        desc_button13.style.button_color = 'lightgreen'
        desc_button14 = Button(description='Base/maximum cell speed', disabled=True, layout=desc_button_layout) 
        desc_button14.style.button_color = 'tan'
        desc_button15 = Button(description='Minimum ECM density required for cell motility', disabled=True, layout=desc_button_layout) 
        desc_button15.style.button_color = 'lightgreen'
        desc_button16 = Button(description='Maximum ECM density allowing cell motility', disabled=True, layout=desc_button_layout) 
        desc_button16.style.button_color = 'tan'
        desc_button17 = Button(description='Ideal ECM density cell motility', disabled=True, layout=desc_button_layout) 
        desc_button17.style.button_color = 'lightgreen'
        desc_button18 = Button(description='Cell directional persistence time', disabled=True, layout=desc_button_layout) 
        desc_button18.style.button_color = 'tan'
        desc_button19 = Button(description='Cell ECM density target', disabled=True, layout=desc_button_layout) 
        desc_button19.style.button_color = 'lightgreen'
        desc_button20 = Button(description='Cell sensitivity to ECM orientation for ECM following type cells', disabled=True, layout=desc_button_layout) 
        desc_button20.style.button_color = 'tan'
        desc_button21 = Button(description="Cell's base sensivity to previous cell direction - cell motility history bias", disabled=True, layout=desc_button_layout) 
        desc_button21.style.button_color = 'lightgreen'
        desc_button22 = Button(description='Chemotaxis bias for non-ECM influenced cell types', disabled=True, layout=desc_button_layout) 
        desc_button22.style.button_color = 'tan'
        desc_button23 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button23.style.button_color = 'lightgreen'
        desc_button24 = Button(description='Chemotactic subsrate decay rate', disabled=True, layout=desc_button_layout) 
        desc_button24.style.button_color = 'tan'
        desc_button25 = Button(description='Allows selection of various cell motility-ECM models', disabled=True, layout=desc_button_layout) 
        desc_button25.style.button_color = 'lightgreen'
        desc_button26 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button26.style.button_color = 'tan'
        desc_button27 = Button(description='Testing and setup parameters', disabled=True, layout=desc_button_layout) 
        desc_button27.style.button_color = 'lightgreen'
        desc_button28 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button28.style.button_color = 'tan'
        desc_button29 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button29.style.button_color = 'lightgreen'
        desc_button30 = Button(description='', disabled=True, layout=desc_button_layout) 
        desc_button30.style.button_color = 'tan'
        desc_button31 = Button(description='Specifies cell, ECM, and domain parameters for consistent unit tests of ECM influenced mechanics and mechanics influence on ECM: 1 for testing set up, 0 for regular set up', disabled=True, layout=desc_button_layout) 
        desc_button31.style.button_color = 'lightgreen'
        desc_button32 = Button(description='Specifies that cells will reset their position to the left boundary and that density is 0.5 for the behavior to ECM march test 1 for march set up, 0 for regular set up', disabled=True, layout=desc_button_layout) 
        desc_button32.style.button_color = 'tan'
        desc_button33 = Button(description='Specifies the initial cell setup: single, random, lesion, cells at y = 0, circle of cells, or cells at left boundary/march', disabled=True, layout=desc_button_layout) 
        desc_button33.style.button_color = 'lightgreen'
        desc_button34 = Button(description='Specifies the initial ECM orietation: random, horizontal, vertical, circular, or starburst', disabled=True, layout=desc_button_layout) 
        desc_button34.style.button_color = 'tan'
        desc_button35 = Button(description='Specifies the chemotatic field gradient orientation: starburst, vertical up, horizontal right, angle, or none', disabled=True, layout=desc_button_layout) 
        desc_button35.style.button_color = 'lightgreen'
        desc_button36 = Button(description='Angle of chemical field gradient orientation, specified in degrees', disabled=True, layout=desc_button_layout) 
        desc_button36.style.button_color = 'tan'

        row1 = [param_name1, self.Model_Initialization, units_button1, desc_button1] 
        row2 = [param_name2, self.ECM_initialization_parameters, units_button2, desc_button2] 
        row3 = [param_name3, self.initial_anisotropy, units_button3, desc_button3] 
        row4 = [param_name4, self.initial_ECM_density, units_button4, desc_button4] 
        row5 = [param_name5, self.Lesion_intialization, units_button5, desc_button5] 
        row6 = [param_name6, self.tumor_radius, units_button6, desc_button6] 
        row7 = [param_name7, self.initial_leader_cell_fraction, units_button7, desc_button7] 
        row8 = [param_name8, self.Cell_properties, units_button8, desc_button8] 
        row9 = [param_name9, self.oxygen_uptake, units_button9, desc_button9] 
        row10 = [param_name10, self.leader_adhesion, units_button10, desc_button10] 
        row11 = [param_name11, self.leader_repulsion, units_button11, desc_button11] 
        row12 = [param_name12, self.follower_adhesion, units_button12, desc_button12] 
        row13 = [param_name13, self.follower_repulsion, units_button13, desc_button13] 
        row14 = [param_name14, self.default_cell_speed, units_button14, desc_button14] 
        row15 = [param_name15, self.rho_L, units_button15, desc_button15] 
        row16 = [param_name16, self.rho_H, units_button16, desc_button16] 
        row17 = [param_name17, self.rho_I, units_button17, desc_button17] 
        row18 = [param_name18, self.default_persistence_time, units_button18, desc_button18] 
        row19 = [param_name19, self.default_ECM_density_target, units_button19, desc_button19] 
        row20 = [param_name20, self.default_ECM_sensitivity, units_button20, desc_button20] 
        row21 = [param_name21, self.default_hysteresis_bias, units_button21, desc_button21] 
        row22 = [param_name22, self.oxygen_migration_bias_for_leaders, units_button22, desc_button22] 
        row23 = [param_name23, self.Microenvironment_parameters, units_button23, desc_button23] 
        row24 = [param_name24, self.chemotactic_substrate_decay_rate, units_button24, desc_button24] 
        row25 = [param_name25, self.Model_Selection_Parameters, units_button25, desc_button25] 
        row26 = [param_name26, self.cell_motility_ECM_interaction_model_selector, units_button26, desc_button26] 
        row27 = [param_name27, self.Testing_and_setup_parameters, units_button27, desc_button27] 
        row28 = [param_name28, self.normalize_ECM_influenced_motility_vector, units_button28, desc_button28] 
        row29 = [param_name29, self.duration_of_uE_conditioning, units_button29, desc_button29] 
        row30 = [param_name30, self.freeze_uE_profile, units_button30, desc_button30] 
        row31 = [param_name31, self.unit_test_setup, units_button31, desc_button31] 
        row32 = [param_name32, self.march_unit_test_setup, units_button32, desc_button32] 
        row33 = [param_name33, self.cell_setup, units_button33, desc_button33] 
        row34 = [param_name34, self.ECM_orientation_setup, units_button34, desc_button34] 
        row35 = [param_name35, self.chemical_field_setup, units_button35, desc_button35] 
        row36 = [param_name36, self.angle_of_chemical_field_gradient, units_button36, desc_button36] 

        box_layout = Layout(display='flex', flex_flow='row', align_items='stretch', width='100%')
        box1 = Box(children=row1, layout=box_layout)
        box2 = Box(children=row2, layout=box_layout)
        box3 = Box(children=row3, layout=box_layout)
        box4 = Box(children=row4, layout=box_layout)
        box5 = Box(children=row5, layout=box_layout)
        box6 = Box(children=row6, layout=box_layout)
        box7 = Box(children=row7, layout=box_layout)
        box8 = Box(children=row8, layout=box_layout)
        box9 = Box(children=row9, layout=box_layout)
        box10 = Box(children=row10, layout=box_layout)
        box11 = Box(children=row11, layout=box_layout)
        box12 = Box(children=row12, layout=box_layout)
        box13 = Box(children=row13, layout=box_layout)
        box14 = Box(children=row14, layout=box_layout)
        box15 = Box(children=row15, layout=box_layout)
        box16 = Box(children=row16, layout=box_layout)
        box17 = Box(children=row17, layout=box_layout)
        box18 = Box(children=row18, layout=box_layout)
        box19 = Box(children=row19, layout=box_layout)
        box20 = Box(children=row20, layout=box_layout)
        box21 = Box(children=row21, layout=box_layout)
        box22 = Box(children=row22, layout=box_layout)
        box23 = Box(children=row23, layout=box_layout)
        box24 = Box(children=row24, layout=box_layout)
        box25 = Box(children=row25, layout=box_layout)
        box26 = Box(children=row26, layout=box_layout)
        box27 = Box(children=row27, layout=box_layout)
        box28 = Box(children=row28, layout=box_layout)
        box29 = Box(children=row29, layout=box_layout)
        box30 = Box(children=row30, layout=box_layout)
        box31 = Box(children=row31, layout=box_layout)
        box32 = Box(children=row32, layout=box_layout)
        box33 = Box(children=row33, layout=box_layout)
        box34 = Box(children=row34, layout=box_layout)
        box35 = Box(children=row35, layout=box_layout)
        box36 = Box(children=row36, layout=box_layout)

        self.tab = VBox([
          box1,
          box2,
          box3,
          box4,
          box5,
          box6,
          box7,
          box8,
          box9,
          box10,
          box11,
          box12,
          box13,
          box14,
          box15,
          box16,
          box17,
          box18,
          box19,
          box20,
          box21,
          box22,
          box23,
          box24,
          box25,
          box26,
          box27,
          box28,
          box29,
          box30,
          box31,
          box32,
          box33,
          box34,
          box35,
          box36,
        ])

    # Populate the GUI widgets with values from the XML
    def fill_gui(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//user_parameters')  # find unique entry point
        self.Model_Initialization.value = float(uep.find('.//Model_Initialization').text)
        self.ECM_initialization_parameters.value = float(uep.find('.//ECM_initialization_parameters').text)
        self.initial_anisotropy.value = float(uep.find('.//initial_anisotropy').text)
        self.initial_ECM_density.value = float(uep.find('.//initial_ECM_density').text)
        self.Lesion_intialization.value = float(uep.find('.//Lesion_intialization').text)
        self.tumor_radius.value = float(uep.find('.//tumor_radius').text)
        self.initial_leader_cell_fraction.value = float(uep.find('.//initial_leader_cell_fraction').text)
        self.Cell_properties.value = float(uep.find('.//Cell_properties').text)
        self.oxygen_uptake.value = float(uep.find('.//oxygen_uptake').text)
        self.leader_adhesion.value = float(uep.find('.//leader_adhesion').text)
        self.leader_repulsion.value = float(uep.find('.//leader_repulsion').text)
        self.follower_adhesion.value = float(uep.find('.//follower_adhesion').text)
        self.follower_repulsion.value = float(uep.find('.//follower_repulsion').text)
        self.default_cell_speed.value = float(uep.find('.//default_cell_speed').text)
        self.rho_L.value = float(uep.find('.//rho_L').text)
        self.rho_H.value = float(uep.find('.//rho_H').text)
        self.rho_I.value = float(uep.find('.//rho_I').text)
        self.default_persistence_time.value = float(uep.find('.//default_persistence_time').text)
        self.default_ECM_density_target.value = float(uep.find('.//default_ECM_density_target').text)
        self.default_ECM_sensitivity.value = float(uep.find('.//default_ECM_sensitivity').text)
        self.default_hysteresis_bias.value = float(uep.find('.//default_hysteresis_bias').text)
        self.oxygen_migration_bias_for_leaders.value = float(uep.find('.//oxygen_migration_bias_for_leaders').text)
        self.Microenvironment_parameters.value = float(uep.find('.//Microenvironment_parameters').text)
        self.chemotactic_substrate_decay_rate.value = float(uep.find('.//chemotactic_substrate_decay_rate').text)
        self.Model_Selection_Parameters.value = float(uep.find('.//Model_Selection_Parameters').text)
        self.cell_motility_ECM_interaction_model_selector.value = (uep.find('.//cell_motility_ECM_interaction_model_selector').text)
        self.Testing_and_setup_parameters.value = float(uep.find('.//Testing_and_setup_parameters').text)
        self.normalize_ECM_influenced_motility_vector.value = ('true' == (uep.find('.//normalize_ECM_influenced_motility_vector').text.lower()) )
        self.duration_of_uE_conditioning.value = float(uep.find('.//duration_of_uE_conditioning').text)
        self.freeze_uE_profile.value = ('true' == (uep.find('.//freeze_uE_profile').text.lower()) )
        self.unit_test_setup.value = int(uep.find('.//unit_test_setup').text)
        self.march_unit_test_setup.value = int(uep.find('.//march_unit_test_setup').text)
        self.cell_setup.value = (uep.find('.//cell_setup').text)
        self.ECM_orientation_setup.value = (uep.find('.//ECM_orientation_setup').text)
        self.chemical_field_setup.value = (uep.find('.//chemical_field_setup').text)
        self.angle_of_chemical_field_gradient.value = float(uep.find('.//angle_of_chemical_field_gradient').text)


    # Read values from the GUI widgets to enable editing XML
    def fill_xml(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//user_parameters')  # find unique entry point
        uep.find('.//Model_Initialization').text = str(self.Model_Initialization.value)
        uep.find('.//ECM_initialization_parameters').text = str(self.ECM_initialization_parameters.value)
        uep.find('.//initial_anisotropy').text = str(self.initial_anisotropy.value)
        uep.find('.//initial_ECM_density').text = str(self.initial_ECM_density.value)
        uep.find('.//Lesion_intialization').text = str(self.Lesion_intialization.value)
        uep.find('.//tumor_radius').text = str(self.tumor_radius.value)
        uep.find('.//initial_leader_cell_fraction').text = str(self.initial_leader_cell_fraction.value)
        uep.find('.//Cell_properties').text = str(self.Cell_properties.value)
        uep.find('.//oxygen_uptake').text = str(self.oxygen_uptake.value)
        uep.find('.//leader_adhesion').text = str(self.leader_adhesion.value)
        uep.find('.//leader_repulsion').text = str(self.leader_repulsion.value)
        uep.find('.//follower_adhesion').text = str(self.follower_adhesion.value)
        uep.find('.//follower_repulsion').text = str(self.follower_repulsion.value)
        uep.find('.//default_cell_speed').text = str(self.default_cell_speed.value)
        uep.find('.//rho_L').text = str(self.rho_L.value)
        uep.find('.//rho_H').text = str(self.rho_H.value)
        uep.find('.//rho_I').text = str(self.rho_I.value)
        uep.find('.//default_persistence_time').text = str(self.default_persistence_time.value)
        uep.find('.//default_ECM_density_target').text = str(self.default_ECM_density_target.value)
        uep.find('.//default_ECM_sensitivity').text = str(self.default_ECM_sensitivity.value)
        uep.find('.//default_hysteresis_bias').text = str(self.default_hysteresis_bias.value)
        uep.find('.//oxygen_migration_bias_for_leaders').text = str(self.oxygen_migration_bias_for_leaders.value)
        uep.find('.//Microenvironment_parameters').text = str(self.Microenvironment_parameters.value)
        uep.find('.//chemotactic_substrate_decay_rate').text = str(self.chemotactic_substrate_decay_rate.value)
        uep.find('.//Model_Selection_Parameters').text = str(self.Model_Selection_Parameters.value)
        uep.find('.//cell_motility_ECM_interaction_model_selector').text = str(self.cell_motility_ECM_interaction_model_selector.value)
        uep.find('.//Testing_and_setup_parameters').text = str(self.Testing_and_setup_parameters.value)
        uep.find('.//normalize_ECM_influenced_motility_vector').text = str(self.normalize_ECM_influenced_motility_vector.value)
        uep.find('.//duration_of_uE_conditioning').text = str(self.duration_of_uE_conditioning.value)
        uep.find('.//freeze_uE_profile').text = str(self.freeze_uE_profile.value)
        uep.find('.//unit_test_setup').text = str(self.unit_test_setup.value)
        uep.find('.//march_unit_test_setup').text = str(self.march_unit_test_setup.value)
        uep.find('.//cell_setup').text = str(self.cell_setup.value)
        uep.find('.//ECM_orientation_setup').text = str(self.ECM_orientation_setup.value)
        uep.find('.//chemical_field_setup').text = str(self.chemical_field_setup.value)
        uep.find('.//angle_of_chemical_field_gradient').text = str(self.angle_of_chemical_field_gradient.value)
