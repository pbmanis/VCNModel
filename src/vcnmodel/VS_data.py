"""
    Vector strength for models with SAM tones, different input configurations.
    17 Aug 2021 version.
    Results are printout from DataTablesVCN after selecting the data runs.
NOTE: This table is automatically written by figures.py and should not be
      directly edited.    pbm
"""

data = """Cell,Configuration,carrierfreq,frequency,dmod,dB,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,AN_VS,AN_phase,AN_phasesd,maxArea,ninputs
2,largestonly,16000.0,200.0,100.0,30.0,0.5408,891,0.5064,1.1087,260.6313,6.4457e-114,0.5326,3.9578,1.1226,174.0100,6
2,removelargest,16000.0,200.0,100.0,30.0,0.8159,10553,5.8412,0.6379,7025.3784,0.0000e+00,0.5326,3.9578,1.1226,174.0100,6
2,largestonly,16000.0,300.0,100.0,30.0,0.5284,780,2.8672,1.1296,217.7634,2.6702e-95,0.5874,5.1563,1.0316,174.0100,6
2,removelargest,16000.0,300.0,100.0,30.0,0.7145,7835,1.9706,0.8200,3999.6088,0.0000e+00,0.5874,5.1563,1.0316,174.0100,6
2,largestonly,16000.0,400.0,100.0,30.0,0.4579,681,5.0935,1.2499,142.7682,9.9213e-63,0.5642,0.0363,1.0699,174.0100,6
2,removelargest,16000.0,400.0,100.0,30.0,0.6430,6687,4.2455,0.9398,2764.4724,0.0000e+00,0.5642,0.0363,1.0699,174.0100,6
2,largestonly,16000.0,500.0,100.0,30.0,0.2932,792,0.7237,1.5665,68.0799,2.7119e-30,0.4906,1.1046,1.1935,174.0100,6
2,removelargest,16000.0,500.0,100.0,30.0,0.5375,5968,0.0552,1.1142,1724.4760,0.0000e+00,0.4906,1.1046,1.1935,174.0100,6
2,largestonly,16000.0,750.0,100.0,30.0,0.1389,902,5.6272,1.9871,17.3956,2.7873e-08,0.2958,3.5093,1.5607,174.0100,6
2,removelargest,16000.0,750.0,100.0,30.0,0.2487,5384,4.9137,1.6683,332.9342,2.5616e-145,0.2958,3.5093,1.5607,174.0100,6
2,largestonly,16000.0,1000.0,100.0,30.0,0.0342,902,3.4488,2.5986,1.0531,3.4885e-01,0.1736,5.6472,1.8712,174.0100,6
2,removelargest,16000.0,1000.0,100.0,30.0,0.0966,5236,3.3098,2.1622,48.8145,6.3113e-22,0.1736,5.6472,1.8712,174.0100,6
2,all,16000.0,200.0,100.0,30.0,0.8670,13899,5.5437,0.5342,10448.9292,0.0000e+00,0.5326,3.9578,1.1226,174.0100,6
2,all,16000.0,300.0,100.0,30.0,0.6941,12478,1.8123,0.8546,6010.9922,0.0000e+00,0.5874,5.1563,1.0316,174.0100,6
2,all,16000.0,400.0,100.0,30.0,0.6360,10680,3.9710,0.9513,4320.5611,0.0000e+00,0.5642,0.0363,1.0699,174.0100,6
2,all,16000.0,500.0,100.0,30.0,0.5247,10169,6.0533,1.1357,2799.8268,0.0000e+00,0.4906,1.1046,1.1935,174.0100,6
2,all,16000.0,750.0,100.0,30.0,0.2353,9459,4.4945,1.7011,523.6641,3.7633e-228,0.2958,3.5093,1.5607,174.0100,6
2,all,16000.0,1000.0,100.0,30.0,0.0872,9223,2.7392,2.2090,70.0809,3.6663e-31,0.1736,5.6472,1.8712,174.0100,6
2,all,16000.0,100.0,100.0,30.0,0.4239,10709,3.0898,1.3101,1924.6627,0.0000e+00,0.3589,2.9892,1.4316,174.0100,6
2,largestonly,16000.0,100.0,100.0,30.0,0.6784,1225,3.7331,0.8810,563.7538,1.4616e-245,0.3589,2.9892,1.4316,174.0100,6
2,removelargest,16000.0,100.0,100.0,30.0,0.6781,8393,3.1840,0.8814,3859.2768,0.0000e+00,0.3589,2.9892,1.4316,174.0100,6
5,largestonly,16000.0,200.0,100.0,30.0,0.6678,153,0.7388,0.8986,68.2304,2.3329e-30,0.5324,3.9581,1.1229,141.0200,7
5,removelargest,16000.0,200.0,100.0,30.0,0.8650,13169,5.6760,0.5387,9852.3430,0.0000e+00,0.5324,3.9581,1.1229,141.0200,7
5,largestonly,16000.0,300.0,100.0,30.0,0.5842,132,3.4819,1.0368,45.0574,2.7027e-20,0.5881,5.1558,1.0304,141.0200,7
5,removelargest,16000.0,300.0,100.0,30.0,0.7139,10395,1.9259,0.8210,5297.5081,0.0000e+00,0.5881,5.1558,1.0304,141.0200,7
5,largestonly,16000.0,400.0,100.0,30.0,0.3875,87,5.5886,1.3770,13.0642,2.1197e-06,0.5637,0.0377,1.0707,141.0200,7
5,removelargest,16000.0,400.0,100.0,30.0,0.6614,8826,4.1548,0.9092,3861.5167,0.0000e+00,0.5637,0.0377,1.0707,141.0200,7
5,largestonly,16000.0,500.0,100.0,30.0,0.1449,111,1.2471,1.9655,2.3314,9.7162e-02,0.4898,1.1056,1.1947,141.0200,7
5,removelargest,16000.0,500.0,100.0,30.0,0.5317,7919,6.2631,1.1239,2239.0351,0.0000e+00,0.4898,1.1056,1.1947,141.0200,7
5,largestonly,16000.0,750.0,100.0,30.0,0.1619,153,5.0961,1.9083,4.0099,1.8135e-02,0.2942,3.5078,1.5642,141.0200,7
5,removelargest,16000.0,750.0,100.0,30.0,0.2547,7046,4.7612,1.6539,457.0450,3.2202e-199,0.2942,3.5078,1.5642,141.0200,7
5,largestonly,16000.0,1000.0,100.0,30.0,0.1382,164,3.2122,1.9896,3.1307,4.3689e-02,0.1734,5.6457,1.8720,141.0200,7
5,removelargest,16000.0,1000.0,100.0,30.0,0.0841,6892,2.8671,2.2252,48.7515,6.7217e-22,0.1734,5.6457,1.8720,141.0200,7
5,all,16000.0,200.0,100.0,30.0,0.9090,14536,5.4487,0.4369,12010.5867,0.0000e+00,0.5324,3.9581,1.1229,141.0200,7
5,all,16000.0,300.0,100.0,30.0,0.7222,13982,1.8053,0.8069,7291.7726,0.0000e+00,0.5881,5.1558,1.0304,141.0200,7
5,all,16000.0,400.0,100.0,30.0,0.6603,11419,3.9723,0.9110,4979.2108,0.0000e+00,0.5637,0.0377,1.0707,141.0200,7
5,all,16000.0,500.0,100.0,30.0,0.5471,10717,6.0633,1.0983,3207.9818,0.0000e+00,0.4898,1.1056,1.1947,141.0200,7
5,all,16000.0,750.0,100.0,30.0,0.2451,9638,4.5087,1.6769,579.0745,3.2446e-252,0.2942,3.5078,1.5642,141.0200,7
5,all,16000.0,1000.0,100.0,30.0,0.0827,9543,2.5947,2.2329,65.2153,4.7570e-29,0.1734,5.6457,1.8720,141.0200,7
5,all,16000.0,100.0,100.0,30.0,0.3759,11208,3.0825,1.3989,1583.6090,0.0000e+00,0.3583,2.9903,1.4328,141.0200,7
5,largestonly,16000.0,100.0,100.0,30.0,0.7605,223,3.9080,0.7399,128.9901,9.5566e-57,0.3583,2.9903,1.4328,141.0200,7
5,removelargest,16000.0,100.0,100.0,30.0,0.5626,9561,3.1011,1.0726,3025.8151,0.0000e+00,0.3583,2.9903,1.4328,141.0200,7
6,all,16000.0,200.0,100.0,30.0,0.9130,14679,5.3940,0.4266,12236.0705,0.0000e+00,0.5326,3.9578,1.1226,159.0100,6
6,all,16000.0,300.0,100.0,30.0,0.7443,15331,1.7196,0.7685,8493.1833,0.0000e+00,0.5874,5.1563,1.0316,159.0100,6
6,all,16000.0,400.0,100.0,30.0,0.6609,11876,3.8504,0.9102,5186.9260,0.0000e+00,0.5642,0.0363,1.0699,159.0100,6
6,all,16000.0,500.0,100.0,30.0,0.5440,11138,5.9266,1.1035,3295.6009,0.0000e+00,0.4906,1.1046,1.1935,159.0100,6
6,all,16000.0,750.0,100.0,30.0,0.2445,10251,4.3073,1.6785,612.7105,8.0022e-267,0.2958,3.5093,1.5607,159.0100,6
6,all,16000.0,1000.0,100.0,30.0,0.0822,10095,2.5266,2.2355,68.1872,2.4359e-30,0.1736,5.6472,1.8712,159.0100,6
6,largestonly,16000.0,200.0,100.0,30.0,0.5710,878,0.3948,1.0587,286.2388,4.8759e-125,0.5326,3.9578,1.1226,159.0100,6
6,removelargest,16000.0,200.0,100.0,30.0,0.8638,13765,5.6101,0.5412,10269.7427,0.0000e+00,0.5326,3.9578,1.1226,159.0100,6
6,largestonly,16000.0,300.0,100.0,30.0,0.5231,762,2.9065,1.1385,208.4842,2.8607e-91,0.5874,5.1563,1.0316,159.0100,6
6,removelargest,16000.0,300.0,100.0,30.0,0.7128,12020,1.8709,0.8229,6106.8316,0.0000e+00,0.5874,5.1563,1.0316,159.0100,6
6,largestonly,16000.0,400.0,100.0,30.0,0.4683,731,5.1691,1.2318,160.3099,2.3894e-70,0.5642,0.0363,1.0699,159.0100,6
6,removelargest,16000.0,400.0,100.0,30.0,0.6488,9966,4.0591,0.9302,4195.0694,0.0000e+00,0.5642,0.0363,1.0699,159.0100,6
6,largestonly,16000.0,500.0,100.0,30.0,0.3116,769,0.9501,1.5272,74.6474,3.8112e-33,0.4906,1.1046,1.1935,159.0100,6
6,removelargest,16000.0,500.0,100.0,30.0,0.5216,9340,6.1597,1.1409,2541.0764,0.0000e+00,0.4906,1.1046,1.1935,159.0100,6
6,largestonly,16000.0,750.0,100.0,30.0,0.1341,843,6.2059,2.0047,15.1499,2.6331e-07,0.2958,3.5093,1.5607,159.0100,6
6,removelargest,16000.0,750.0,100.0,30.0,0.2370,8478,4.6700,1.6968,476.3420,1.3407e-207,0.2958,3.5093,1.5607,159.0100,6
6,largestonly,16000.0,1000.0,100.0,30.0,0.1081,918,5.5600,2.1093,10.7285,2.1911e-05,0.1736,5.6472,1.8712,159.0100,6
6,removelargest,16000.0,1000.0,100.0,30.0,0.0858,8329,2.7378,2.2161,61.3261,2.3250e-27,0.1736,5.6472,1.8712,159.0100,6
6,all,16000.0,100.0,100.0,30.0,0.3394,11644,3.1395,1.4701,1341.1787,0.0000e+00,0.3589,2.9892,1.4316,159.0100,6
6,largestonly,16000.0,100.0,100.0,30.0,0.6568,1177,3.7606,0.9169,507.7443,3.0862e-221,0.3589,2.9892,1.4316,159.0100,6
6,removelargest,16000.0,100.0,100.0,30.0,0.4546,10419,3.1194,1.2557,2152.9703,0.0000e+00,0.3589,2.9892,1.4316,159.0100,6
9,all,16000.0,200.0,100.0,30.0,0.8643,15194,5.2807,0.5401,11349.4645,0.0000e+00,0.5329,3.9575,1.1219,245.7600,8
9,all,16000.0,300.0,100.0,30.0,0.7244,18802,1.6448,0.8030,9866.8626,0.0000e+00,0.5896,5.1522,1.0280,245.7600,8
9,all,16000.0,400.0,100.0,30.0,0.5703,16314,3.9265,1.0597,5306.8081,0.0000e+00,0.5631,0.0401,1.0716,245.7600,8
9,all,16000.0,500.0,100.0,30.0,0.4642,15419,5.8748,1.2390,3322.0553,0.0000e+00,0.4900,1.1059,1.1944,245.7600,8
9,all,16000.0,750.0,100.0,30.0,0.1998,15404,4.3167,1.7946,615.0947,7.3754e-268,0.2954,3.5085,1.5617,245.7600,8
9,all,16000.0,1000.0,100.0,30.0,0.0718,15289,2.3133,2.2953,78.7709,6.1690e-35,0.1735,5.6418,1.8718,245.7600,8
9,largestonly,16000.0,200.0,100.0,30.0,0.6689,10844,6.0119,0.8968,4852.0294,0.0000e+00,0.5329,3.9575,1.1219,245.7600,8
9,removelargest,16000.0,200.0,100.0,30.0,0.8901,14514,5.5810,0.4825,11499.3723,0.0000e+00,0.5329,3.9575,1.1219,245.7600,8
9,largestonly,16000.0,300.0,100.0,30.0,0.5373,9913,2.0598,1.1146,2861.7280,0.0000e+00,0.5896,5.1522,1.0280,245.7600,8
9,removelargest,16000.0,300.0,100.0,30.0,0.6791,14053,1.9914,0.8798,6480.6144,0.0000e+00,0.5896,5.1522,1.0280,245.7600,8
9,largestonly,16000.0,400.0,100.0,30.0,0.5299,10059,4.2064,1.1270,2824.3867,0.0000e+00,0.5631,0.0401,1.0716,245.7600,8
9,removelargest,16000.0,400.0,100.0,30.0,0.5989,12306,4.2045,1.0126,4413.8320,0.0000e+00,0.5631,0.0401,1.0716,245.7600,8
9,largestonly,16000.0,500.0,100.0,30.0,0.4135,10123,0.0962,1.3290,1730.9135,0.0000e+00,0.4900,1.1059,1.1944,245.7600,8
9,removelargest,16000.0,500.0,100.0,30.0,0.4789,12014,0.0659,1.2135,2755.4008,0.0000e+00,0.4900,1.1059,1.1944,245.7600,8
9,largestonly,16000.0,750.0,100.0,30.0,0.2290,10103,4.9477,1.7169,530.0007,6.6624e-231,0.2954,3.5085,1.5617,245.7600,8
9,removelargest,16000.0,750.0,100.0,30.0,0.2034,11412,4.8929,1.7846,472.2398,8.1073e-206,0.2954,3.5085,1.5617,245.7600,8
9,largestonly,16000.0,1000.0,100.0,30.0,0.1162,10087,3.3094,2.0750,136.0840,7.9336e-60,0.1735,5.6418,1.8718,245.7600,8
9,removelargest,16000.0,1000.0,100.0,30.0,0.0702,11438,3.0909,2.3051,56.3401,3.4026e-25,0.1735,5.6418,1.8718,245.7600,8
9,all,16000.0,100.0,100.0,30.0,0.2192,14043,3.5927,1.7422,674.9008,7.8397e-294,0.3581,2.9907,1.4330,245.7600,8
9,largestonly,16000.0,100.0,100.0,30.0,0.3685,9738,3.4760,1.4131,1322.1089,0.0000e+00,0.3581,2.9907,1.4330,245.7600,8
9,removelargest,16000.0,100.0,100.0,30.0,0.3152,11782,3.1696,1.5197,1170.2557,0.0000e+00,0.3581,2.9907,1.4330,245.7600,8
10,all,16000.0,200.0,100.0,30.0,0.9322,14704,5.5132,0.3748,12776.7072,0.0000e+00,0.5331,3.9554,1.1217,121.5700,10
10,all,16000.0,300.0,100.0,30.0,0.7561,13806,1.9766,0.7477,7893.0608,0.0000e+00,0.5890,5.1523,1.0290,121.5700,10
10,all,16000.0,400.0,100.0,30.0,0.6899,10804,4.2250,0.8617,5142.1548,0.0000e+00,0.5634,0.0407,1.0712,121.5700,10
10,all,16000.0,500.0,100.0,30.0,0.5700,9611,0.1387,1.0603,3122.3094,0.0000e+00,0.4897,1.1048,1.1949,121.5700,10
10,all,16000.0,750.0,100.0,30.0,0.2576,7970,5.0123,1.6471,528.8014,2.2105e-230,0.2957,3.5089,1.5610,121.5700,10
10,all,16000.0,1000.0,100.0,30.0,0.0713,7849,3.1999,2.2979,39.9554,4.4422e-18,0.1751,5.6422,1.8667,121.5700,10
10,largestonly,16000.0,200.0,100.0,30.0,0.7577,23,1.0581,0.7449,13.2055,1.8405e-06,0.5331,3.9554,1.1217,121.5700,10
10,removelargest,16000.0,200.0,100.0,30.0,0.8975,13771,5.7113,0.4650,11093.4308,0.0000e+00,0.5331,3.9554,1.1217,121.5700,10
10,largestonly,16000.0,300.0,100.0,30.0,0.7372,17,3.6229,0.7808,9.2399,9.7084e-05,0.5890,5.1523,1.0290,121.5700,10
10,removelargest,16000.0,300.0,100.0,30.0,0.7555,10530,2.0797,0.7488,6010.8696,0.0000e+00,0.5890,5.1523,1.0290,121.5700,10
10,largestonly,16000.0,400.0,100.0,30.0,0.3966,17,0.5037,1.3600,2.6740,6.8974e-02,0.5634,0.0407,1.0712,121.5700,10
10,removelargest,16000.0,400.0,100.0,30.0,0.6901,8498,4.4208,0.8614,4046.5400,0.0000e+00,0.5634,0.0407,1.0712,121.5700,10
10,largestonly,16000.0,500.0,100.0,30.0,0.2556,10,3.7705,1.6518,0.6531,5.2043e-01,0.4897,1.1048,1.1949,121.5700,10
10,removelargest,16000.0,500.0,100.0,30.0,0.5835,7087,0.3488,1.0380,2413.0401,0.0000e+00,0.4897,1.1048,1.1949,121.5700,10
10,largestonly,16000.0,750.0,100.0,30.0,0.2770,23,6.0261,1.6023,1.7652,1.7116e-01,0.2957,3.5089,1.5610,121.5700,10
10,removelargest,16000.0,750.0,100.0,30.0,0.2570,5708,5.3427,1.6483,377.1439,1.6162e-164,0.2957,3.5089,1.5610,121.5700,10
10,largestonly,16000.0,1000.0,100.0,30.0,0.4019,24,4.8170,1.3503,3.8762,2.0730e-02,0.1751,5.6422,1.8667,121.5700,10
10,removelargest,16000.0,1000.0,100.0,30.0,0.0741,5486,3.6838,2.2811,30.1627,7.9529e-14,0.1751,5.6422,1.8667,121.5700,10
10,all,16000.0,100.0,100.0,30.0,0.4399,10626,3.0590,1.2816,2056.0745,0.0000e+00,0.3583,2.9885,1.4327,121.5700,10
10,largestonly,16000.0,100.0,100.0,30.0,0.8267,45,4.0006,0.6169,30.7568,4.3902e-14,0.3583,2.9885,1.4327,121.5700,10
10,removelargest,16000.0,100.0,100.0,30.0,0.6029,9405,3.1055,1.0060,3418.6044,0.0000e+00,0.3583,2.9885,1.4327,121.5700,10
11,all,16000.0,200.0,100.0,30.0,0.8359,15382,5.1474,0.5987,10748.4356,0.0000e+00,0.5324,3.9581,1.1229,213.0900,7
11,all,16000.0,300.0,100.0,30.0,0.7379,19569,1.4451,0.7797,10654.0328,0.0000e+00,0.5881,5.1558,1.0304,213.0900,7
11,all,16000.0,400.0,100.0,30.0,0.5690,17180,3.7228,1.0620,5562.0579,0.0000e+00,0.5637,0.0377,1.0707,213.0900,7
11,all,16000.0,500.0,100.0,30.0,0.4553,16096,5.6397,1.2543,3337.4007,0.0000e+00,0.4898,1.1056,1.1947,213.0900,7
11,all,16000.0,750.0,100.0,30.0,0.2081,16305,4.0181,1.7719,705.8752,2.7688e-307,0.2942,3.5078,1.5642,213.0900,7
11,all,16000.0,1000.0,100.0,30.0,0.0734,16173,1.9155,2.2855,87.1581,1.4052e-38,0.1734,5.6457,1.8720,213.0900,7
11,largestonly,16000.0,200.0,100.0,30.0,0.6444,10481,5.9123,0.9375,4352.2892,0.0000e+00,0.5324,3.9581,1.1229,213.0900,7
11,removelargest,16000.0,200.0,100.0,30.0,0.8680,14602,5.3880,0.5321,11001.8997,0.0000e+00,0.5324,3.9581,1.1229,213.0900,7
11,largestonly,16000.0,300.0,100.0,30.0,0.5495,9671,1.8865,1.0942,2920.5575,0.0000e+00,0.5881,5.1558,1.0304,213.0900,7
11,removelargest,16000.0,300.0,100.0,30.0,0.6599,15678,1.7422,0.9117,6827.9009,0.0000e+00,0.5881,5.1558,1.0304,213.0900,7
11,largestonly,16000.0,400.0,100.0,30.0,0.5244,9857,4.0327,1.1362,2710.9333,0.0000e+00,0.5637,0.0377,1.0707,213.0900,7
11,removelargest,16000.0,400.0,100.0,30.0,0.5636,13835,3.8672,1.0709,4394.9453,0.0000e+00,0.5637,0.0377,1.0707,213.0900,7
11,largestonly,16000.0,500.0,100.0,30.0,0.4255,9839,6.0746,1.3073,1781.1896,0.0000e+00,0.4898,1.1056,1.1947,213.0900,7
11,removelargest,16000.0,500.0,100.0,30.0,0.4610,13530,5.9061,1.2444,2875.8129,0.0000e+00,0.4898,1.1056,1.1947,213.0900,7
11,largestonly,16000.0,750.0,100.0,30.0,0.2204,9846,4.5798,1.7392,478.1308,2.2409e-208,0.2942,3.5078,1.5642,213.0900,7
11,removelargest,16000.0,750.0,100.0,30.0,0.2099,13289,4.3064,1.7670,585.4747,5.3896e-255,0.2942,3.5078,1.5642,213.0900,7
11,largestonly,16000.0,1000.0,100.0,30.0,0.1175,9844,2.8255,2.0693,136.0003,8.6264e-60,0.1734,5.6457,1.8720,213.0900,7
11,removelargest,16000.0,1000.0,100.0,30.0,0.0656,13220,2.2826,2.3343,56.8656,2.0119e-25,0.1734,5.6457,1.8720,213.0900,7
11,all,16000.0,100.0,100.0,30.0,0.2096,14481,3.6492,1.7678,636.2862,4.6175e-277,0.3583,2.9903,1.4328,213.0900,7
11,largestonly,16000.0,100.0,100.0,30.0,0.3736,9524,3.4293,1.4033,1329.1126,0.0000e+00,0.3583,2.9903,1.4328,213.0900,7
11,removelargest,16000.0,100.0,100.0,30.0,0.2429,12814,3.2688,1.6822,756.2832,0.0000e+00,0.3583,2.9903,1.4328,213.0900,7
13,all,16000.0,200.0,100.0,30.0,0.9175,14830,5.2948,0.4149,12484.7572,0.0000e+00,0.5324,3.9581,1.1229,146.5700,7
13,all,16000.0,300.0,100.0,30.0,0.7235,16788,1.6647,0.8046,8787.7184,0.0000e+00,0.5881,5.1558,1.0304,146.5700,7
13,all,16000.0,400.0,100.0,30.0,0.6189,13902,3.8602,0.9796,5324.8252,0.0000e+00,0.5637,0.0377,1.0707,146.5700,7
13,all,16000.0,500.0,100.0,30.0,0.5218,13312,5.9161,1.1406,3624.2261,0.0000e+00,0.4898,1.1056,1.1947,146.5700,7
13,all,16000.0,750.0,100.0,30.0,0.2364,12681,4.3416,1.6984,708.6872,1.6637e-308,0.2942,3.5078,1.5642,146.5700,7
13,all,16000.0,1000.0,100.0,30.0,0.0799,12570,2.3540,2.2481,80.2629,1.3876e-35,0.1734,5.6457,1.8720,146.5700,7
13,largestonly,16000.0,200.0,100.0,30.0,0.5357,4148,0.0977,1.1173,1190.3009,0.0000e+00,0.5324,3.9581,1.1229,146.5700,7
13,removelargest,16000.0,200.0,100.0,30.0,0.8763,14127,5.5305,0.5139,10848.2502,0.0000e+00,0.5324,3.9581,1.1229,146.5700,7
13,largestonly,16000.0,300.0,100.0,30.0,0.5444,4049,2.4013,1.1029,1199.8026,0.0000e+00,0.5881,5.1558,1.0304,146.5700,7
13,removelargest,16000.0,300.0,100.0,30.0,0.6909,12809,1.8594,0.8600,6113.9475,0.0000e+00,0.5881,5.1558,1.0304,146.5700,7
13,largestonly,16000.0,400.0,100.0,30.0,0.4872,4076,4.6669,1.1993,967.3236,0.0000e+00,0.5637,0.0377,1.0707,146.5700,7
13,removelargest,16000.0,400.0,100.0,30.0,0.6336,11060,4.0416,0.9554,4439.7056,0.0000e+00,0.5637,0.0377,1.0707,146.5700,7
13,largestonly,16000.0,500.0,100.0,30.0,0.4045,4158,0.5189,1.3454,680.3649,3.3210e-296,0.4898,1.1056,1.1947,146.5700,7
13,removelargest,16000.0,500.0,100.0,30.0,0.5116,10528,6.1822,1.1578,2755.3703,0.0000e+00,0.4898,1.1056,1.1947,146.5700,7
13,largestonly,16000.0,750.0,100.0,30.0,0.1765,4193,5.7067,1.8625,130.6368,1.8414e-57,0.2942,3.5078,1.5642,146.5700,7
13,removelargest,16000.0,750.0,100.0,30.0,0.2461,9661,4.7067,1.6745,585.2660,6.6407e-255,0.2942,3.5078,1.5642,146.5700,7
13,largestonly,16000.0,1000.0,100.0,30.0,0.0754,4155,4.4238,2.2735,23.6475,5.3705e-11,0.1734,5.6457,1.8720,146.5700,7
13,removelargest,16000.0,1000.0,100.0,30.0,0.0731,9518,2.8311,2.2872,50.8917,7.9072e-23,0.1734,5.6457,1.8720,146.5700,7
13,all,16000.0,100.0,100.0,30.0,0.2566,12648,3.1931,1.6495,832.5603,0.0000e+00,0.3583,2.9903,1.4328,146.5700,7
13,largestonly,16000.0,100.0,100.0,30.0,0.5952,4563,3.6084,1.0186,1616.6431,0.0000e+00,0.3583,2.9903,1.4328,146.5700,7
13,removelargest,16000.0,100.0,100.0,30.0,0.4128,10817,3.0665,1.3302,1843.7028,0.0000e+00,0.3583,2.9903,1.4328,146.5700,7
17,all,16000.0,200.0,100.0,30.0,0.8169,15055,5.1443,0.6359,10047.5649,0.0000e+00,0.5324,3.9581,1.1229,278.3200,7
17,all,16000.0,300.0,100.0,30.0,0.7479,18777,1.2561,0.7622,10503.8615,0.0000e+00,0.5881,5.1558,1.0304,278.3200,7
17,all,16000.0,400.0,100.0,30.0,0.5989,16078,3.4787,1.0125,5767.7353,0.0000e+00,0.5637,0.0377,1.0707,278.3200,7
17,all,16000.0,500.0,100.0,30.0,0.4633,14817,5.2972,1.2405,3180.3891,0.0000e+00,0.4898,1.1056,1.1947,278.3200,7
17,all,16000.0,750.0,100.0,30.0,0.2388,15347,3.5355,1.6924,875.1763,0.0000e+00,0.2942,3.5078,1.5642,278.3200,7
17,all,16000.0,1000.0,100.0,30.0,0.1110,15233,1.3446,2.0967,187.7449,2.9068e-82,0.1734,5.6457,1.8720,278.3200,7
17,largestonly,16000.0,200.0,100.0,30.0,0.6804,11492,5.6378,0.8775,5320.5611,0.0000e+00,0.5324,3.9581,1.1229,278.3200,7
17,removelargest,16000.0,200.0,100.0,30.0,0.7903,13817,5.4056,0.6860,8630.7662,0.0000e+00,0.5324,3.9581,1.1229,278.3200,7
17,largestonly,16000.0,300.0,100.0,30.0,0.5884,11445,1.7165,1.0300,3961.7575,0.0000e+00,0.5881,5.1558,1.0304,278.3200,7
17,removelargest,16000.0,300.0,100.0,30.0,0.6688,14546,1.5327,0.8970,6505.9914,0.0000e+00,0.5881,5.1558,1.0304,278.3200,7
17,largestonly,16000.0,400.0,100.0,30.0,0.5035,10853,3.7164,1.1715,2751.4934,0.0000e+00,0.5637,0.0377,1.0707,278.3200,7
17,removelargest,16000.0,400.0,100.0,30.0,0.5567,12976,3.6681,1.0823,4021.9473,0.0000e+00,0.5637,0.0377,1.0707,278.3200,7
17,largestonly,16000.0,500.0,100.0,30.0,0.4688,11055,5.6578,1.2309,2429.9360,0.0000e+00,0.4898,1.1056,1.1947,278.3200,7
17,removelargest,16000.0,500.0,100.0,30.0,0.4644,12531,5.5517,1.2385,2702.7341,0.0000e+00,0.4898,1.1056,1.1947,278.3200,7
17,largestonly,16000.0,750.0,100.0,30.0,0.2463,11042,4.0135,1.6740,669.8287,1.2505e-291,0.2942,3.5078,1.5642,278.3200,7
17,removelargest,16000.0,750.0,100.0,30.0,0.2302,12535,3.8609,1.7138,664.4852,2.6167e-289,0.2942,3.5078,1.5642,278.3200,7
17,largestonly,16000.0,1000.0,100.0,30.0,0.1282,11078,2.1489,2.0268,182.1583,7.7565e-80,0.1734,5.6457,1.8720,278.3200,7
17,removelargest,16000.0,1000.0,100.0,30.0,0.1016,12527,1.7179,2.1384,129.4068,6.3002e-57,0.1734,5.6457,1.8720,278.3200,7
17,all,16000.0,100.0,100.0,30.0,0.2317,14033,3.5062,1.7102,753.1649,0.0000e+00,0.3583,2.9903,1.4328,278.3200,7
17,largestonly,16000.0,100.0,100.0,30.0,0.3279,10532,3.5235,1.4933,1132.5565,0.0000e+00,0.3583,2.9903,1.4328,278.3200,7
17,removelargest,16000.0,100.0,100.0,30.0,0.3158,12078,3.3360,1.5184,1204.2663,0.0000e+00,0.3583,2.9903,1.4328,278.3200,7
18,all,16000.0,200.0,100.0,30.0,0.8892,14725,5.2685,0.4847,11641.4769,0.0000e+00,0.5329,3.9575,1.1219,273.2200,8
18,all,16000.0,300.0,100.0,30.0,0.7456,17064,1.4902,0.7662,9487.1214,0.0000e+00,0.5896,5.1522,1.0280,273.2200,8
18,all,16000.0,400.0,100.0,30.0,0.6301,13643,3.6504,0.9612,5416.2312,0.0000e+00,0.5631,0.0401,1.0716,273.2200,8
18,all,16000.0,500.0,100.0,30.0,0.5143,12614,5.5876,1.1531,3337.0954,0.0000e+00,0.4900,1.1059,1.1944,273.2200,8
18,all,16000.0,750.0,100.0,30.0,0.2437,12291,3.8575,1.6803,730.1742,7.7513e-318,0.2954,3.5085,1.5617,273.2200,8
18,all,16000.0,1000.0,100.0,30.0,0.0993,12300,1.7246,2.1491,121.3447,1.9982e-53,0.1735,5.6418,1.8718,273.2200,8
18,largestonly,16000.0,200.0,100.0,30.0,0.6653,10660,5.8985,0.9028,4718.0340,0.0000e+00,0.5329,3.9575,1.1219,273.2200,8
18,removelargest,16000.0,200.0,100.0,30.0,0.8723,13373,5.5878,0.5227,10175.9517,0.0000e+00,0.5329,3.9575,1.1219,273.2200,8
18,largestonly,16000.0,300.0,100.0,30.0,0.5462,9669,1.9250,1.0998,2884.5967,0.0000e+00,0.5896,5.1522,1.0280,273.2200,8
18,removelargest,16000.0,300.0,100.0,30.0,0.7312,11196,1.8207,0.7912,5986.4324,0.0000e+00,0.5896,5.1522,1.0280,273.2200,8
18,largestonly,16000.0,400.0,100.0,30.0,0.5380,9808,4.0333,1.1135,2838.4795,0.0000e+00,0.5631,0.0401,1.0716,273.2200,8
18,removelargest,16000.0,400.0,100.0,30.0,0.6540,9082,4.0327,0.9216,3884.3394,0.0000e+00,0.5631,0.0401,1.0716,273.2200,8
18,largestonly,16000.0,500.0,100.0,30.0,0.4274,9805,6.1113,1.3039,1790.7471,0.0000e+00,0.4900,1.1059,1.1944,273.2200,8
18,removelargest,16000.0,500.0,100.0,30.0,0.5477,8183,6.0856,1.0973,2454.8536,0.0000e+00,0.4900,1.1059,1.1944,273.2200,8
18,largestonly,16000.0,750.0,100.0,30.0,0.2357,9835,4.5978,1.7000,546.5919,4.1509e-238,0.2954,3.5085,1.5617,273.2200,8
18,removelargest,16000.0,750.0,100.0,30.0,0.2571,7307,4.5861,1.6481,483.1444,1.4896e-210,0.2954,3.5085,1.5617,273.2200,8
18,largestonly,16000.0,1000.0,100.0,30.0,0.1290,9795,2.8726,2.0240,162.8898,1.8108e-71,0.1735,5.6418,1.8718,273.2200,8
18,removelargest,16000.0,1000.0,100.0,30.0,0.0954,7229,2.6152,2.1677,65.8107,2.6229e-29,0.1735,5.6418,1.8718,273.2200,8
18,all,16000.0,100.0,100.0,30.0,0.2899,12542,3.2777,1.5737,1053.8417,0.0000e+00,0.3581,2.9907,1.4330,273.2200,8
18,largestonly,16000.0,100.0,100.0,30.0,0.3867,9567,3.4063,1.3784,1430.8933,0.0000e+00,0.3581,2.9907,1.4330,273.2200,8
18,removelargest,16000.0,100.0,100.0,30.0,0.5174,9915,3.1076,1.1480,2654.4536,0.0000e+00,0.3581,2.9907,1.4330,273.2200,8
30,all,16000.0,200.0,100.0,30.0,0.9162,14631,5.5278,0.4185,12280.3316,0.0000e+00,0.5328,3.9559,1.1222,172.1400,9
30,all,16000.0,300.0,100.0,30.0,0.7122,14189,2.0057,0.8238,7197.6780,0.0000e+00,0.5891,5.1531,1.0288,172.1400,9
30,all,16000.0,400.0,100.0,30.0,0.6293,11749,4.2265,0.9625,4652.6149,0.0000e+00,0.5628,0.0398,1.0721,172.1400,9
30,all,16000.0,500.0,100.0,30.0,0.5167,11178,0.1337,1.1491,2984.7522,0.0000e+00,0.4898,1.1046,1.1949,172.1400,9
30,all,16000.0,750.0,100.0,30.0,0.2269,10209,4.9947,1.7225,525.4144,6.5378e-229,0.2952,3.5111,1.5621,172.1400,9
30,all,16000.0,1000.0,100.0,30.0,0.0666,9981,3.3154,2.3274,44.3272,5.6096e-20,0.1746,5.6425,1.8682,172.1400,9
30,largestonly,16000.0,200.0,100.0,30.0,0.5578,1586,0.5758,1.0805,493.5075,4.7035e-215,0.5328,3.9559,1.1222,172.1400,9
30,removelargest,16000.0,200.0,100.0,30.0,0.8777,13161,5.8390,0.5107,10139.3746,0.0000e+00,0.5328,3.9559,1.1222,172.1400,9
30,largestonly,16000.0,300.0,100.0,30.0,0.5343,1483,2.9551,1.1197,423.3348,1.4061e-184,0.5891,5.1531,1.0288,172.1400,9
30,removelargest,16000.0,300.0,100.0,30.0,0.7322,9487,2.1917,0.7896,5085.9286,0.0000e+00,0.5891,5.1531,1.0288,172.1400,9
30,largestonly,16000.0,400.0,100.0,30.0,0.4756,1415,5.2472,1.2191,320.1262,9.3535e-140,0.5628,0.0398,1.0721,172.1400,9
30,removelargest,16000.0,400.0,100.0,30.0,0.6610,7954,4.5902,0.9099,3475.3138,0.0000e+00,0.5628,0.0398,1.0721,172.1400,9
30,largestonly,16000.0,500.0,100.0,30.0,0.3364,1533,1.2254,1.4761,173.5132,4.4074e-76,0.4898,1.1046,1.1949,172.1400,9
30,removelargest,16000.0,500.0,100.0,30.0,0.5545,6726,0.5082,1.0860,2067.8803,0.0000e+00,0.4898,1.1046,1.1949,172.1400,9
30,largestonly,16000.0,750.0,100.0,30.0,0.1139,1650,0.3560,2.0845,21.3970,5.0982e-10,0.2952,3.5111,1.5621,172.1400,9
30,removelargest,16000.0,750.0,100.0,30.0,0.2407,5604,5.6345,1.6877,324.6656,9.9892e-142,0.2952,3.5111,1.5621,172.1400,9
30,largestonly,16000.0,1000.0,100.0,30.0,0.0357,1689,5.2489,2.5816,2.1543,1.1598e-01,0.1746,5.6425,1.8682,172.1400,9
30,removelargest,16000.0,1000.0,100.0,30.0,0.0636,5443,4.1516,2.3471,22.0490,2.6560e-10,0.1746,5.6425,1.8682,172.1400,9
30,all,16000.0,100.0,100.0,30.0,0.3383,11557,3.1142,1.4724,1322.4198,0.0000e+00,0.3584,2.9892,1.4326,172.1400,9
30,largestonly,16000.0,100.0,100.0,30.0,0.6655,1950,3.7507,0.9024,863.7281,0.0000e+00,0.3584,2.9892,1.4326,172.1400,9
30,removelargest,16000.0,100.0,100.0,30.0,0.6260,9184,3.1403,0.9679,3599.1101,0.0000e+00,0.3584,2.9892,1.4326,172.1400,9
"""
