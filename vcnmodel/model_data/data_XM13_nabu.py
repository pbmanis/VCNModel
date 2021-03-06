"""
Data structures for model_run
XM13 with nabu  channels

"""
#     nabu _vshift   4.3    [1]     4.3    [1]     4.3    [1]    4.3    [1]    4.3    [1]

ChannelData = u"""

    This table describes the REFERENCE ion channel densities (and voltage shifts if necessary)
    for different cell types based on the Xie and Manis 2013 models for mouse.

    The REFERENCE values are applied to "point" models, and to the soma of
    compartmental models.
    The names of the mechanisms must match a channel mechanism (Neuron .mod files)
    and the following _(gbar, vshift, etc) must match an attribute of that channel
    that can be accessed.

    -----------------------------------------------------------------------------------------------------------------------------------
                   II             II-I           I-c           I-II          I-t       
                                                                           
    nabu_gbar      1000.  [1]     1000.  [1]     3000.  [1]    1000.  [2]    3000.  [1] 
    kht_gbar       58.0   [3]     58.0   [1]     500.0  [1]    150.0  [2]    500.0  [1] 
    klt_gbar       80.0   [1]     14.0   [1]     0.0    [1]    20.0   [2]    0.0    [1] 
    ka_gbar        0.0    [1]     0.0    [1]     0.0    [1]    0.0    [2]    125.0  [1] 
    ihvcn_gbar     30.0   [1]     30.0   [1]     18.0   [1]    2.0    [2]    18.0   [1] 
    leak_gbar      2.0    [1]     2.0    [1]     8.0    [1]    2.0    [2]    8.0    [1] 
    leak_erev      -65    [1]     -65    [1]     -65    [1]    -65    [2]    -65    [1] 
    na_type        nabu   [1]     nabu   [1]     nabu   [1]    nabu   [1]    nabu   [1] 
    ih_type        ihvcn  [1]     ihvcn  [1]     ihvcn  [1]    ihvcn  [2]    ihvcn  [1] 
    soma_Cap       26.0   [1]     26.0   [1]     25.0   [1]    26.0   [2]    25.0   [1] 
    e_k            -84.   [1]     -84.   [1]     -84.   [1]    -84.   [2]    -84.   [1] 
    e_na           50.    [1]     50.    [1]     50.    [1]    50.    [2]    50.    [1] 
    ih_eh          -43.   [1]     -43.   [1]     -43.   [1]    -43.   [2]    -43.   [1]
    nabu_vshift    -00.           0.             0.            0.            0.

    -----------------------------------------------------------------------------------------------------------------------------------

    [1] Uses channels from Rothman and Manis, 2003
        Conductances are for Mouse bushy cells
        Xie and Manis, 2013
        Age "adult", Temperature=34C
        Units are nS.
        nabu_vshift: for stability... 

    [2] Rothman and Manis, 2003, model I-II
        Some low-voltage K current, based on observations of
        a single spike near threshold and regular firing for higher
        currents (Xie and Manis, 2017)
    
    [3] Upped delayed rectifier to get faster repolarizaition



    """

#    nabu _vshift    4.3  [2]   4.3 [2]              0.0 [2]            4.3 [2]           4.3 [2]     4.3 [2]     0.0  [2]         0.0  [2]            0.0 [2]

ChannelCompartments = u"""

    Densities for decoration of a cell
    This table describes the ion channel densities **relative** to somatic densities,
    e.g., relative to REFERENCE densities in the table XM13_channels (ChannelData).
    and voltage shifts, for different compartments of the specified neuron,
    Conductances will be calculated from the Model derived from Xie and Manis 2013 for mouse
    (data table: XM13_channels).

    NOTE: unmyelinatedaxon and initialsegment are equivalent in George's models, but only "unmyelinatedaxon" is actually used.
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------
                   axon       unmyelinatedaxon     myelinatedaxon     initialsegment    hillock     soma        dendrite         primarydendrite    secondarydendrite

    nabu_gbar      3.0 [1]    15. [1]              0.0 [1]            0.0 [1]           9.0 [1]     2.5 [1]     0.5 [1]          0.25 [1]           0.25 [1]
    kht_gbar       1.0 [1]    1.0 [1]              0.01 [1]           2.0 [1]           1.0 [1]     1.0 [1]     1.0 [1]          0.5 [1]            0.25 [1]
    klt_gbar       1.0 [1]    1.0 [1]              0.01 [1]           1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.25 [1]
    ihvcn_gbar     0.0 [1]    0.0 [1]              0.0 [1]            0.5 [1]           0.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]
    leak_gbar      1.0 [1]    0.25 [1]             0.25e-3 [1]        1.0 [1]           1.0 [1]     1.0 [1]     0.5 [1]          0.5 [1]            0.5 [1]
    leak_erev      -65. [1]   -65. [1]             -65. [1]           -65. [1]          -65. [1]    -65. [1]    -65. [1]         -65. [1]           -65. [1]
    na_type        nabu       nabu                 nabu               nabu              nabu        nabu        nabu             nabu               nabu  
    ih_type        ihvcn      ihvcn                ihvcn              ihvcn             ihvcn       ihvcn       ihvcn            ihvcn              ihvcn
    nabu_vshift    0          0                    0                  0                 0           0           0                0                  0
    -------------------------------------------------------------------------------------------------------------------------------------------------------------------

    [1] Scaling is relative to soma scaling. Numbers are estimates based on general distribution from literature on cortical neurons.
    [2] Set to 0 (was 4.3 in original model). Matches original Barela et al (2006) scaling.
    [3] for the shift, the value is added to the soma value. 

    """
