"""
    Vector strength for models with SAM tones, different input configurations.
    17 Aug 2021 version.
    Results are a printout from DataTablesVCN after selecting the data runs.
Run started at: 2023-04-04-14.12.06
WARNING: This table is automatically written by figures.py generate_VS_data_file
       and should not be directly edited.
To Regenerate:
   After running the simulationns, select the runs to be added.
   Use the 'Print File Info' button for each cells
    Copy the text in the 'Reports' dock in DataTablesVCN
  into a 'VS_datasets_xxdB.py' file, where xx is the sound pressure level.
Then select 'VS-SAMTone-no figure' in DataTables, and 'Create Figure.
  No figure will be generated, but vector strength will be calculated
   and the VS_data_xxdB.py file will be created.
The VS_data_xxdB.py file holds all of the vector-strength information, in a text format,
   and is read by the plotting programs.
--pbm 2014-2022
"""

data = """Cell,Filename,Configuration,carrierfreq,frequency,dmod,dB,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,VS_mean,VS_SD,VS_Ns,VS_groups,AN_VS,AN_phase,AN_phasesd,SAC_AN,SAC_Bu,SAC_AN_HW,SAC_Bu_HW,maxArea,ninputs
2,runANPSTH-all-2021-12-02.13-30-16,all,16000.0,0050.0,100.0,15.0,0.6834,5979,3.2273,0.8726,2792.0486,0.0000e+00,0.6835,0.0075,5979,10,0.7144,3.0978,0.8201,1.1423,4.0205,0.008307,0.002158,174.0100,6
2,runANPSTH-all-2021-12-03.01-55-35,all,16000.0,0100.0,100.0,15.0,0.9492,7545,4.0503,0.3229,6797.7356,0.0000e+00,0.9492,0.0055,7545,10,0.7853,3.6814,0.6952,1.3233,7.0690,0.003539,0.001255,174.0100,6
2,runANPSTH-all-2021-12-03.18-47-11,all,16000.0,0200.0,100.0,15.0,0.8978,14113,0.4500,0.4643,11375.7401,0.0000e+00,0.8979,0.0041,14113,10,0.8124,4.9521,0.6447,1.3809,4.1223,0.001709,0.001091,174.0100,6
2,runANPSTH-all-2021-12-06.04-07-38,all,16000.0,0300.0,100.0,15.0,0.8160,12344,2.8762,0.6377,8219.6944,0.0000e+00,0.8160,0.0073,12344,10,0.8004,6.1556,0.6673,1.3145,2.8367,0.001212,0.001085,174.0100,6
2,runANPSTH-all-2021-12-06.16-41-20,all,16000.0,0400.0,100.0,15.0,0.7552,10207,4.9445,0.7494,5820.9479,0.0000e+00,0.7553,0.0103,10207,10,0.7712,0.9697,0.7209,1.2166,2.4623,0.000990,0.000923,174.0100,6
2,runANPSTH-all-2021-12-07.05-04-38,all,16000.0,0500.0,100.0,15.0,0.6676,9213,0.6586,0.8989,4106.3920,0.0000e+00,0.6679,0.0123,9213,10,0.7302,1.9711,0.7931,1.1301,2.0170,0.000860,0.000851,174.0100,6
2,runANPSTH-all-2021-12-07.22-18-16,all,16000.0,0750.0,100.0,15.0,0.3910,7707,5.2680,1.3705,1178.0632,0.0000e+00,0.3915,0.0231,7707,10,0.5557,4.1830,1.0839,1.0383,1.3593,0.000639,0.000628,174.0100,6
2,runANPSTH-all-2021-12-08.10-37-14,all,16000.0,1000.0,100.0,15.0,0.1333,7423,3.3005,2.0074,131.9699,4.8551e-58,0.1373,0.0210,7423,10,0.3085,6.2002,1.5337,1.0449,1.0826,0.000477,0.000644,174.0100,6
2,runANPSTH-largestonly-2021-12-02.13-50-32,largestonly,16000.0,0050.0,100.0,15.0,0.7756,1074,3.3095,0.7129,646.0865,2.5596e-281,0.7775,0.0331,1074,10,0.7144,3.0978,0.8201,1.1400,3.8343,0.008328,0.003893,174.0100,6
2,runANPSTH-largestonly-2021-12-03.02-15-55,largestonly,16000.0,0100.0,100.0,15.0,0.8661,1890,4.7007,0.5362,1417.7014,0.0000e+00,0.8672,0.0132,1890,10,0.7853,3.6814,0.6952,1.3237,3.6460,0.003535,0.002588,174.0100,6
2,runANPSTH-largestonly-2021-12-05.16-10-32,largestonly,16000.0,0200.0,100.0,15.0,0.7950,1114,1.2696,0.6773,704.1288,1.5877e-306,0.7957,0.0240,1114,10,0.8124,4.9521,0.6447,1.3835,2.9904,0.001705,0.001770,174.0100,6
2,runANPSTH-largestonly-2021-12-06.04-28-20,largestonly,16000.0,0300.0,100.0,15.0,0.7668,900,3.5540,0.7287,529.2189,1.4559e-230,0.7679,0.0237,900,10,0.8004,6.1556,0.6673,1.3166,2.8108,0.001212,0.001446,174.0100,6
2,runANPSTH-largestonly-2021-12-06.17-03-21,largestonly,16000.0,0400.0,100.0,15.0,0.6992,900,5.8548,0.8460,439.9472,8.5777e-192,0.7011,0.0451,900,10,0.7712,0.9697,0.7209,1.2171,2.4782,0.000990,0.001118,174.0100,6
2,runANPSTH-largestonly-2021-12-07.05-25-14,largestonly,16000.0,0500.0,100.0,15.0,0.5808,932,1.7859,1.0424,314.3946,2.8852e-137,0.5809,0.0602,932,10,0.7302,1.9711,0.7931,1.1318,2.1289,0.000858,0.000874,174.0100,6
2,runANPSTH-largestonly-2021-12-07.22-39-22,largestonly,16000.0,0750.0,100.0,15.0,0.3234,1139,0.2594,1.5026,119.1185,1.8514e-52,0.3264,0.0371,1139,10,0.5557,4.1830,1.0839,1.0375,1.4951,0.000640,0.000768,174.0100,6
2,runANPSTH-largestonly-2021-12-08.10-59-01,largestonly,16000.0,1000.0,100.0,15.0,0.1385,1065,5.3589,1.9885,20.4223,1.3512e-09,0.1561,0.0463,1065,10,0.3085,6.2002,1.5337,1.0426,1.3289,0.000477,0.000485,174.0100,6
2,runANPSTH-removelargest-2021-12-02.14-11-28,removelargest,16000.0,0050.0,100.0,15.0,0.7555,4526,3.1426,0.7489,2583.0766,0.0000e+00,0.7557,0.0099,4526,10,0.7144,3.0978,0.8201,1.1393,4.7040,0.008320,0.002469,174.0100,6
2,runANPSTH-removelargest-2021-12-03.02-36-23,removelargest,16000.0,0100.0,100.0,15.0,0.9455,7242,4.2269,0.3349,6473.5186,0.0000e+00,0.9455,0.0029,7242,10,0.7853,3.6814,0.6952,1.3230,5.9987,0.003539,0.001471,174.0100,6
2,runANPSTH-removelargest-2021-12-05.16-30-55,removelargest,16000.0,0200.0,100.0,15.0,0.8761,11129,0.7640,0.5144,8541.4647,0.0000e+00,0.8761,0.0044,11129,10,0.8124,4.9521,0.6447,1.3844,3.5579,0.001706,0.001304,174.0100,6
2,runANPSTH-removelargest-2021-12-06.04-48-52,removelargest,16000.0,0300.0,100.0,15.0,0.8236,7869,3.0423,0.6231,5337.0907,0.0000e+00,0.8238,0.0076,7869,10,0.8004,6.1556,0.6673,1.3137,2.9867,0.001213,0.001042,174.0100,6
2,runANPSTH-removelargest-2021-12-06.17-24-30,removelargest,16000.0,0400.0,100.0,15.0,0.7621,6613,5.1985,0.7371,3840.7407,0.0000e+00,0.7627,0.0078,6613,10,0.7712,0.9697,0.7209,1.2179,2.4506,0.000989,0.000940,174.0100,6
2,runANPSTH-removelargest-2021-12-07.05-45-20,removelargest,16000.0,0500.0,100.0,15.0,0.6805,5451,0.9276,0.8775,2524.0072,0.0000e+00,0.6810,0.0105,5451,10,0.7302,1.9711,0.7931,1.1309,2.0841,0.000859,0.000830,174.0100,6
2,runANPSTH-removelargest-2021-12-07.23-00-10,removelargest,16000.0,0750.0,100.0,15.0,0.4068,3961,5.6408,1.3412,655.5604,1.9666e-285,0.4073,0.0331,3961,10,0.5557,4.1830,1.0839,1.0427,1.4129,0.000636,0.000638,174.0100,6
2,runANPSTH-removelargest-2021-12-08.11-19-41,removelargest,16000.0,1000.0,100.0,15.0,0.1472,3766,3.7221,1.9577,81.5546,3.8132e-36,0.1520,0.0329,3766,10,0.3085,6.2002,1.5337,1.0440,1.1614,0.000476,0.000501,174.0100,6