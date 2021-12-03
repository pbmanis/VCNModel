"""
Data structures for model_run
XM13 with nacn  channels

"""
# add_table_data('XM13_channels_nacncoop', row_key='field', col_key='model_type',
#                species='mouse', data=

ChannelData = u"""

This table describes the REFERENCE ion channel densities (and voltage shifts if necessary)
for different cell types based on the Xie and Manis 2013 models for mouse, but using
the nacncoop mechanism (coooperative sodium channels)

!!!!!!!!!!!! USAGE OF THIS TABLE SHOULD BE CONSIDERED EXPERIMENTAL !!!!!!!!!!!!!!

The REFERENCE values are applied to "point" models, and to the soma of
compartmental models.
The names of the mechanisms must match a channel mechanism (Neuron .mod files)
and the following _(gbar, vshift, etc) must match an attribute of that channel
that can be accessed.

-----------------------------------------------------------------------------------------------------------------------------------
                 II             II-I           I-c           I-II          I-t       
                                                                                     
nacncoop_gbar    2500.  [4]     1000.  [4]     1000.  [4]    1000.  [4]    1000.  [4] 
kht_gbar         58.0   [1]     58.0   [1]     500.0  [1]    150.0  [2]    500.0  [1] 
klt_gbar         80.0   [1]     20.0   [1]     0.0    [1]    14.0   [3]    0.0    [1] 
ka_gbar          0.0    [1]     0.0    [1]     0.0    [1]    0.0    [2]    125.0  [1] 
ihvcn_gbar       30.0   [1]     30.0   [1]     18.0   [1]    2.0    [2]    18.0   [1] 
leak_gbar        2.0    [1]     2.0    [1]     8.0    [1]    2.0    [2]    8.0    [1] 
leak_erev        -65    [1]     -65    [1]     -65    [1]    -65    [2]    -65    [1] 
na_type          nacncoop [1]   nacncoop  [1]  nacncoop [1]  nacncoop [3]  nacncoop [1] 
ih_type          ihvcn  [1]     ihvcn  [1]     ihvcn  [1]    ihvcn  [2]    ihvcn  [1] 
soma_Cap         26.0   [1]     26.0   [1]     25.0   [1]    25.0   [2]    25.0   [1] 
nacncoop_vshift  0.0    [1]     0.0    [1]     0.0    [1]    0.0    [1]    0.0    [1]
e_k              -84    [1]     -84    [1]     -84    [1]    -70    [3]    -84    [1] 
e_na             50.    [1]     50.    [1]     50.    [1]    55.    [3]    50.    [1] 
ih_eh            -43    [1]     -43    [1]     -43    [1]    -43    [2]    -43    [1] 

-----------------------------------------------------------------------------------------------------------------------------------

[1] Uses channels from Xie and Manis, 2013
    Age "adult", Temperature=34C
    Units are nS.

[2] Rothman and Manis, 2003, model I-II
    Some low-voltage K current, based on observations of
    a single spike near threshold and regular firing for higher
    currents (Xie and Manis, 2017)

[3] These values for the I-II (dstellate) are from the original checkpoint test
    for cnmodel 12/2017. 

[4] nav11 channels were used in original Xie and Manis (2013) ms, 
    However, this version uses cooperative na channels for faster activation

"""

# add_table_data('XM13_channels_nacncoop_compartments', row_key='parameter', col_key='compartment',
#                species='mouse', model_type='II', data=u"""

ChannelCompartments=u"""
!!!!!!!!!!!! USAGE OF THIS TABLE SHOULD BE CONSIDERED EXPERIMENTAL !!!!!!!!!!!!!!

This table describes the ion channel densities relative to somatic densities,
e.g., relative to REFERENCE densities in the table XM13_nacncoop_channels.
and voltage shifts, for different compartments of the specified neuron,
Conductances will be calculated from the Model derived from Xie and Manis 2013 for mouse
Note:
4/10/2019
reduced dendrite Na to 0 (was 0.5). 
------------------------------------------------------------------------------------------------------------------------------------------------------------------
                   axon           unmyelinatedaxon     myelinatedaxon     initialsegment    hillock     soma        dendrite         primarydendrite    secondarydendrite
                                                                                                                                                                                                      
nacncoop_gbar      1.0 [1]        5.0 [1]              1.0 [1]            5.0 [1]           5.0 [1]     1.0 [1]     0.0 [1]          0.0 [1]            0.0 [1]       
kht_gbar           1.0 [1]        2.0 [1]              0.01 [1]           2.0 [1]           2.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
klt_gbar           1.0 [1]        1.0 [1]              0.01 [1]           1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]       
ihvcn_gbar         0.0 [1]        0.0 [1]              0.0 [1]            0.5 [1]           0.0 [1]     1.0 [1]     0.5[1]           0.5 [1]            0.5 [1]       
leak_gbar          1.0 [1]        0.25 [1]             0.25e-3 [1]        1.0 [1]           1.0 [1]     1.0 [1]     5.0 [1]          0.5 [1]            0.5 [1]       
leak_erev          -65. [1]       -65. [1]             -65. [1]           -65. [1]          -65. [1]    -65. [1]    -65. [1]         -65. [1]           -65. [1]      
nacncoop_vshift    0.0  [1]       0.0  [1]             0.0 [1]            0.0  [1]          0.0  [1]    0.0 [1]     0.0 [1]          0.0 [1]            0.0 [1]       
na_type            nacncoop       nacncoop             nacncoop           nacncoop          nacncoop    nacncoop    nacncoop            nacncoop              nacncoop
ih_type            ihvcn          ihvcn                ihvcn              ihvcn             ihvcn       ihvcn       ihvcn            ihvcn              ihvcn                            
--------------------------------------------------------------------------------------------------------------------------------------------------------------------

[1] Scaling is relative to soma scaling. Numbers are estimates based on general distribution from literature on cortical neurons.


"""

# ***** END OF XM13_Channels for nacncoop version of model