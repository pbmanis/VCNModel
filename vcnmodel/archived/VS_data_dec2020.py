"""
    Vector strength for models with SAM tones, different input configurations.
    9 Dec 2020 version.
    Results are printout from DataTablesVCN after selecting the data runs.
    pbm
"""

data = """Cell,Configuration,frequency,dmod,dB,VectorStrength,SpikeCount,phase,phasesd,Rayleigh,RayleighP,AN_VS,maxArea,ninputs
2,all,1000.0,100.0,20.0,0.1341,7684,0.5047,0.3190,138.1787,9.7665e-61,0.2927,174.0100,6
2,largestonly,1000.0,100.0,20.0,0.1006,797,0.8669,0.3411,8.0603,3.1584e-04,0.2927,174.0100,6
2,removelargest,1000.0,100.0,20.0,0.1242,4029,0.5884,0.3251,62.1233,1.0476e-27,0.2927,174.0100,6
2,all,750.0,100.0,20.0,0.3883,8009,0.8162,0.2189,1207.8518,0.0000e+00,0.5288,174.0100,6
2,largestonly,750.0,100.0,20.0,0.2432,859,0.0669,0.2676,50.8130,8.5547e-23,0.5288,174.0100,6
2,removelargest,750.0,100.0,20.0,0.3756,4309,0.8906,0.2227,608.0147,8.7618e-265,0.5288,174.0100,6
2,all,500.0,100.0,20.0,0.6682,9756,0.0840,0.1429,4356.0487,0.0000e+00,0.7061,174.0100,6
2,largestonly,500.0,100.0,20.0,0.5570,545,0.2700,0.1722,169.0568,3.7982e-74,0.7061,174.0100,6
2,removelargest,500.0,100.0,20.0,0.6798,5754,0.1322,0.1398,2658.7190,0.0000e+00,0.7061,174.0100,6
2,all,400.0,100.0,20.0,0.7546,10900,0.7529,0.1194,6206.6167,0.0000e+00,0.7509,174.0100,6
2,largestonly,400.0,100.0,20.0,0.6311,533,0.9260,0.1527,212.3147,6.2071e-93,0.7509,174.0100,6
2,removelargest,400.0,100.0,20.0,0.7611,6991,0.7955,0.1176,4049.9209,0.0000e+00,0.7509,174.0100,6
2,all,300.0,100.0,20.0,0.8152,14115,0.4120,0.1018,9379.0624,0.0000e+00,0.7860,174.0100,6
2,largestonly,300.0,100.0,20.0,0.7475,577,0.5513,0.1214,322.4348,9.2971e-141,0.7860,174.0100,6
2,removelargest,300.0,100.0,20.0,0.8147,8919,0.4418,0.1019,5920.1142,0.0000e+00,0.7860,174.0100,6
2,all,200.0,100.0,20.0,0.9193,14670,0.9873,0.0653,12398.9859,0.0000e+00,0.7981,174.0100,6
2,largestonly,200.0,100.0,20.0,0.8123,848,0.1626,0.1026,559.5179,1.0103e-243,0.7981,174.0100,6
2,removelargest,200.0,100.0,20.0,0.8834,12812,0.0438,0.0792,9999.4674,0.0000e+00,0.7981,174.0100,6
2,all,100.0,100.0,20.0,0.8880,7908,0.5802,0.0776,6236.3968,0.0000e+00,0.7054,174.0100,6
2,largestonly,100.0,100.0,20.0,0.8861,1527,0.6805,0.0783,1199.0719,0.0000e+00,0.7054,174.0100,6
2,removelargest,100.0,100.0,20.0,0.9530,7408,0.6009,0.0494,6727.6078,0.0000e+00,0.7054,174.0100,6

5,all,1000.0,100.0,20.0,0.1429,69,0.4977,0.3140,1.4090,2.4439e-01,0.2925,141.0200,7
5,largestonly,1000.0,100.0,20.0,nan,0,nan,nan,nan,1.0000e+00,0.2925,141.0200,7
5,removelargest,1000.0,100.0,20.0,0.3317,12,0.9880,0.2365,1.3199,2.6715e-01,0.2925,141.0200,7
5,all,750.0,100.0,20.0,0.6634,141,0.7596,0.1442,62.0602,1.1158e-27,0.5285,141.0200,7
5,largestonly,750.0,100.0,20.0,nan,0,nan,nan,nan,1.0000e+00,0.5285,141.0200,7
5,removelargest,750.0,100.0,20.0,0.6530,23,0.7771,0.1469,9.8079,5.5018e-05,0.5285,141.0200,7
5,all,500.0,100.0,20.0,0.8679,417,0.0404,0.0847,314.0724,3.9821e-137,0.7067,141.0200,7
5,largestonly,500.0,100.0,20.0,nan,0,nan,nan,nan,1.0000e+00,0.7067,141.0200,7
5,removelargest,500.0,100.0,20.0,0.8747,84,0.0713,0.0823,64.2732,1.2204e-28,0.7067,141.0200,7
5,all,400.0,100.0,20.0,0.9046,702,0.7174,0.0713,574.4619,3.2687e-250,0.7503,141.0200,7
5,largestonly,400.0,100.0,20.0,nan,0,nan,nan,nan,1.0000e+00,0.7503,141.0200,7
5,removelargest,400.0,100.0,20.0,0.9147,134,0.7568,0.0672,112.1064,2.0549e-49,0.7503,141.0200,7
5,all,300.0,100.0,20.0,0.9314,2176,0.3982,0.0600,1887.6898,0.0000e+00,0.7864,141.0200,7
5,largestonly,300.0,100.0,20.0,nan,0,nan,nan,nan,1.0000e+00,0.7864,141.0200,7
5,removelargest,300.0,100.0,20.0,0.9280,349,0.4037,0.0615,300.5283,3.0355e-131,0.7864,141.0200,7
5,all,200.0,100.0,20.0,0.9561,9569,0.0464,0.0477,8746.9614,0.0000e+00,0.7982,141.0200,7
5,largestonly,200.0,100.0,20.0,nan,0,nan,nan,nan,1.0000e+00,0.7982,141.0200,7
5,removelargest,200.0,100.0,20.0,0.9593,3578,0.0589,0.0459,3292.5126,0.0000e+00,0.7982,141.0200,7
5,all,100.0,100.0,20.0,0.9875,5138,0.6070,0.0252,5010.6937,0.0000e+00,0.7050,141.0200,7
5,largestonly,100.0,100.0,20.0,nan,0,nan,nan,nan,1.0000e+00,0.7050,141.0200,7
5,removelargest,100.0,100.0,20.0,0.9878,2433,0.6164,0.0250,2373.8390,0.0000e+00,0.7050,141.0200,7

6,all,1000.0,100.0,20.0,0.1255,10758,0.4853,0.3242,169.5205,2.3889e-74,0.2927,159.0100,6
6,largestonly,1000.0,100.0,20.0,0.0998,1753,0.8496,0.3417,17.4436,2.6568e-08,0.2927,159.0100,6
6,removelargest,1000.0,100.0,20.0,0.1222,8970,0.5660,0.3263,133.9311,6.8306e-59,0.2927,159.0100,6
6,all,750.0,100.0,20.0,0.3982,10937,0.7937,0.2160,1734.2127,0.0000e+00,0.5288,159.0100,6
6,largestonly,750.0,100.0,20.0,0.2767,1779,0.1144,0.2551,136.2105,6.9906e-60,0.5288,159.0100,6
6,removelargest,750.0,100.0,20.0,0.3789,9205,0.8428,0.2217,1321.7176,0.0000e+00,0.5288,159.0100,6
6,all,500.0,100.0,20.0,0.6964,12308,0.0524,0.1354,5969.0860,0.0000e+00,0.7061,159.0100,6
6,largestonly,500.0,100.0,20.0,0.5537,1471,0.2708,0.1731,450.9630,1.4102e-196,0.7061,159.0100,6
6,removelargest,500.0,100.0,20.0,0.6770,10807,0.0915,0.1406,4953.5518,0.0000e+00,0.7061,159.0100,6
6,all,400.0,100.0,20.0,0.7804,13283,0.7287,0.1121,8089.3181,0.0000e+00,0.7509,159.0100,6
6,largestonly,400.0,100.0,20.0,0.6640,1384,0.9162,0.1440,610.1620,1.0233e-265,0.7509,159.0100,6
6,removelargest,400.0,100.0,20.0,0.7562,11437,0.7606,0.1190,6540.4586,0.0000e+00,0.7509,159.0100,6
6,all,300.0,100.0,20.0,0.8366,18179,0.3926,0.0951,12724.4952,0.0000e+00,0.7860,159.0100,6
6,largestonly,300.0,100.0,20.0,0.7467,1391,0.5426,0.1217,775.4962,0.0000e+00,0.7860,159.0100,6
6,removelargest,300.0,100.0,20.0,0.8121,14537,0.4198,0.1027,9587.4597,0.0000e+00,0.7860,159.0100,6
6,all,200.0,100.0,20.0,0.9435,14948,0.9616,0.0543,13306.1879,0.0000e+00,0.7981,159.0100,6
6,largestonly,200.0,100.0,20.0,0.7819,1854,0.1530,0.1117,1133.3755,0.0000e+00,0.7981,159.0100,6
6,removelargest,200.0,100.0,20.0,0.9145,14707,0.9972,0.0673,12299.3569,0.0000e+00,0.7981,159.0100,6
6,all,100.0,100.0,20.0,0.8013,8408,0.5755,0.1059,5398.0752,0.0000e+00,0.7054,159.0100,6
6,largestonly,100.0,100.0,20.0,0.8810,2717,0.6795,0.0801,2109.0036,0.0000e+00,0.7054,159.0100,6
6,removelargest,100.0,100.0,20.0,0.8759,7970,0.5848,0.0819,6114.7339,0.0000e+00,0.7054,159.0100,6

9,all,1000.0,100.0,20.0,0.1233,15522,0.4400,0.3256,235.9310,3.4392e-103,0.2913,245.7600,8
9,largestonly,1000.0,100.0,20.0,0.1962,9970,0.5505,0.2873,383.6754,2.3543e-167,0.2913,245.7600,8
9,removelargest,1000.0,100.0,20.0,0.1179,11949,0.5583,0.3291,166.1258,7.1201e-73,0.2913,245.7600,8
9,all,750.0,100.0,20.0,0.3392,15756,0.7777,0.2340,1812.5191,0.0000e+00,0.5280,245.7600,8
9,largestonly,750.0,100.0,20.0,0.4334,9986,0.8441,0.2058,1875.8990,0.0000e+00,0.5280,245.7600,8
9,removelargest,750.0,100.0,20.0,0.3308,11966,0.8657,0.2367,1309.2827,0.0000e+00,0.5280,245.7600,8
9,all,500.0,100.0,20.0,0.6185,15755,0.0412,0.1560,6026.0696,0.0000e+00,0.7066,245.7600,8
9,largestonly,500.0,100.0,20.0,0.6220,10475,0.1016,0.1551,4053.1474,0.0000e+00,0.7066,245.7600,8
9,removelargest,500.0,100.0,20.0,0.6342,13148,0.1090,0.1519,5288.6924,0.0000e+00,0.7066,245.7600,8
9,all,400.0,100.0,20.0,0.7196,17835,0.7333,0.1291,9234.7983,0.0000e+00,0.7500,245.7600,8
9,largestonly,400.0,100.0,20.0,0.7038,10265,0.7629,0.1334,5085.2308,0.0000e+00,0.7500,245.7600,8
9,removelargest,400.0,100.0,20.0,0.7253,13839,0.7754,0.1276,7280.1329,0.0000e+00,0.7500,245.7600,8
9,all,300.0,100.0,20.0,0.8454,20834,0.3650,0.0922,14889.5863,0.0000e+00,0.7870,245.7600,8
9,largestonly,300.0,100.0,20.0,0.7184,10488,0.4386,0.1294,5412.4498,0.0000e+00,0.7870,245.7600,8
9,removelargest,300.0,100.0,20.0,0.8119,16766,0.4291,0.1027,11052.9692,0.0000e+00,0.7870,245.7600,8
9,all,200.0,100.0,20.0,0.9563,14998,0.9422,0.0476,13714.9082,0.0000e+00,0.7986,245.7600,8
9,largestonly,200.0,100.0,20.0,0.8403,11718,0.0496,0.0939,8274.6708,0.0000e+00,0.7986,245.7600,8
9,removelargest,200.0,100.0,20.0,0.9401,14922,0.9867,0.0559,13188.3818,0.0000e+00,0.7986,245.7600,8
9,all,100.0,100.0,20.0,0.5664,10431,0.6058,0.1697,3345.9266,0.0000e+00,0.7050,245.7600,8
9,largestonly,100.0,100.0,20.0,0.7392,8067,0.6352,0.1237,4407.4824,0.0000e+00,0.7050,245.7600,8
9,removelargest,100.0,100.0,20.0,0.7677,8591,0.5864,0.1157,5063.6998,0.0000e+00,0.7050,245.7600,8

10,all,1000.0,100.0,20.0,0.1152,8033,0.6448,0.3309,106.5673,5.2291e-47,0.2904,121.5700,10
10,largestonly,1000.0,100.0,20.0,0.3763,38,0.8249,0.2225,5.3818,4.5995e-03,0.2904,121.5700,10
10,removelargest,1000.0,100.0,20.0,0.1207,5666,0.6895,0.3273,82.5948,1.3475e-36,0.2904,121.5700,10
10,all,750.0,100.0,20.0,0.3727,8396,0.9231,0.2236,1166.2694,0.0000e+00,0.5268,121.5700,10
10,largestonly,750.0,100.0,20.0,0.4467,33,0.0041,0.2021,6.5845,1.3816e-03,0.5268,121.5700,10
10,removelargest,750.0,100.0,20.0,0.3797,6119,0.9661,0.2215,882.3069,0.0000e+00,0.5268,121.5700,10
10,all,500.0,100.0,20.0,0.6878,11172,0.1422,0.1377,5285.8619,0.0000e+00,0.7063,121.5700,10
10,largestonly,500.0,100.0,20.0,0.1772,10,0.3072,0.2961,0.3138,7.3064e-01,0.7063,121.5700,10
10,removelargest,500.0,100.0,20.0,0.6839,8702,0.1750,0.1387,4069.8994,0.0000e+00,0.7063,121.5700,10
10,all,400.0,100.0,20.0,0.7755,12412,0.7964,0.1135,7464.5752,0.0000e+00,0.7505,121.5700,10
10,largestonly,400.0,100.0,20.0,0.8635,8,0.1017,0.0862,5.9654,2.5659e-03,0.7505,121.5700,10
10,removelargest,400.0,100.0,20.0,0.7710,10336,0.8264,0.1148,6143.4172,0.0000e+00,0.7505,121.5700,10
10,all,300.0,100.0,20.0,0.8388,16466,0.4424,0.0944,11584.1207,0.0000e+00,0.7871,121.5700,10
10,largestonly,300.0,100.0,20.0,0.8875,29,0.6853,0.0777,22.8437,1.1997e-10,0.7871,121.5700,10
10,removelargest,300.0,100.0,20.0,0.8219,13054,0.4642,0.0997,8818.2689,0.0000e+00,0.7871,121.5700,10
10,all,200.0,100.0,20.0,0.9571,14982,0.9853,0.0471,13725.4288,0.0000e+00,0.7989,121.5700,10
10,largestonly,200.0,100.0,20.0,0.8409,60,0.2470,0.0937,42.4224,3.7687e-19,0.7989,121.5700,10
10,removelargest,200.0,100.0,20.0,0.9363,14795,0.0188,0.0577,12970.0842,0.0000e+00,0.7989,121.5700,10
10,all,100.0,100.0,20.0,0.9051,7869,0.5773,0.0711,6446.6562,0.0000e+00,0.7054,121.5700,10
10,largestonly,100.0,100.0,20.0,0.9325,139,0.7369,0.0595,120.8625,3.2364e-53,0.7054,121.5700,10
10,removelargest,100.0,100.0,20.0,0.9563,7627,0.5876,0.0476,6975.6391,0.0000e+00,0.7054,121.5700,10

11,all,1000.0,100.0,20.0,0.1337,16197,0.3972,0.3193,289.3739,2.1208e-126,0.2925,213.0900,7
11,largestonly,1000.0,100.0,20.0,0.2031,9683,0.5084,0.2842,399.3985,3.4949e-174,0.2925,213.0900,7
11,removelargest,1000.0,100.0,20.0,0.1154,13390,0.4694,0.3307,178.3147,3.6218e-78,0.2925,213.0900,7
11,all,750.0,100.0,20.0,0.3436,16469,0.7354,0.2326,1944.4460,0.0000e+00,0.5285,213.0900,7
11,largestonly,750.0,100.0,20.0,0.4225,9732,0.8141,0.2089,1736.9107,0.0000e+00,0.5285,213.0900,7
11,removelargest,750.0,100.0,20.0,0.3369,13538,0.7899,0.2348,1536.9615,0.0000e+00,0.5285,213.0900,7
11,all,500.0,100.0,20.0,0.6211,16185,0.0094,0.1553,6244.3126,0.0000e+00,0.7067,213.0900,7
11,largestonly,500.0,100.0,20.0,0.6239,10131,0.0772,0.1546,3944.0759,0.0000e+00,0.7067,213.0900,7
11,removelargest,500.0,100.0,20.0,0.6215,14120,0.0528,0.1552,5454.4162,0.0000e+00,0.7067,213.0900,7
11,all,400.0,100.0,20.0,0.7147,18635,0.7107,0.1304,9519.3430,0.0000e+00,0.7503,213.0900,7
11,largestonly,400.0,100.0,20.0,0.7012,10060,0.7517,0.1341,4946.3936,0.0000e+00,0.7503,213.0900,7
11,removelargest,400.0,100.0,20.0,0.7066,15089,0.7377,0.1327,7532.9119,0.0000e+00,0.7503,213.0900,7
11,all,300.0,100.0,20.0,0.8508,21072,0.3384,0.0905,15251.4255,0.0000e+00,0.7864,213.0900,7
11,largestonly,300.0,100.0,20.0,0.7289,10095,0.4192,0.1266,5363.7716,0.0000e+00,0.7864,213.0900,7
11,removelargest,300.0,100.0,20.0,0.8073,17871,0.3922,0.1041,11647.7934,0.0000e+00,0.7864,213.0900,7
11,all,200.0,100.0,20.0,0.9537,14996,0.9265,0.0490,13639.7690,0.0000e+00,0.7982,213.0900,7
11,largestonly,200.0,100.0,20.0,0.8375,11406,0.0436,0.0948,7999.7457,0.0000e+00,0.7982,213.0900,7
11,removelargest,200.0,100.0,20.0,0.9307,14888,0.9645,0.0603,12897.3162,0.0000e+00,0.7982,213.0900,7
11,all,100.0,100.0,20.0,0.5202,10881,0.6073,0.1820,2944.2895,0.0000e+00,0.7050,213.0900,7
11,largestonly,100.0,100.0,20.0,0.7577,7897,0.6281,0.1186,4533.2596,0.0000e+00,0.7050,213.0900,7
11,removelargest,100.0,100.0,20.0,0.6519,9445,0.5929,0.1472,4014.1570,0.0000e+00,0.7050,213.0900,7

13,all,1000.0,100.0,20.0,0.1286,12396,0.4902,0.3224,204.9281,1.0020e-89,0.2925,146.5700,7
13,largestonly,1000.0,100.0,20.0,0.1405,4728,0.7693,0.3153,93.2891,3.0552e-41,0.2925,146.5700,7
13,removelargest,1000.0,100.0,20.0,0.1305,9499,0.5553,0.3212,161.8739,5.0009e-71,0.2925,146.5700,7
13,all,750.0,100.0,20.0,0.3740,12727,0.7965,0.2232,1780.6406,0.0000e+00,0.5285,146.5700,7
13,largestonly,750.0,100.0,20.0,0.3218,4722,0.0101,0.2397,488.9214,4.6147e-213,0.5285,146.5700,7
13,removelargest,750.0,100.0,20.0,0.3751,9799,0.8536,0.2229,1378.6190,0.0000e+00,0.5285,146.5700,7
13,all,500.0,100.0,20.0,0.6663,13940,0.0570,0.1434,6188.8946,0.0000e+00,0.7067,146.5700,7
13,largestonly,500.0,100.0,20.0,0.5809,4591,0.2042,0.1659,1549.1275,0.0000e+00,0.7067,146.5700,7
13,removelargest,500.0,100.0,20.0,0.6550,11499,0.0991,0.1464,4933.6618,0.0000e+00,0.7067,146.5700,7
13,all,400.0,100.0,20.0,0.7450,14921,0.7346,0.1221,8280.7142,0.0000e+00,0.7503,146.5700,7
13,largestonly,400.0,100.0,20.0,0.6498,4414,0.8621,0.1478,1863.5939,0.0000e+00,0.7503,146.5700,7
13,removelargest,400.0,100.0,20.0,0.7405,12099,0.7646,0.1234,6633.9842,0.0000e+00,0.7503,146.5700,7
13,all,300.0,100.0,20.0,0.8363,18822,0.3846,0.0952,13165.3402,0.0000e+00,0.7864,146.5700,7
13,largestonly,300.0,100.0,20.0,0.7469,4406,0.4960,0.1216,2458.0644,0.0000e+00,0.7864,146.5700,7
13,removelargest,300.0,100.0,20.0,0.8074,14681,0.4211,0.1041,9570.6638,0.0000e+00,0.7864,146.5700,7
13,all,200.0,100.0,20.0,0.9505,14979,0.9512,0.0507,13533.4179,0.0000e+00,0.7982,146.5700,7
13,largestonly,200.0,100.0,20.0,0.7805,4678,0.1154,0.1120,2849.7310,0.0000e+00,0.7982,146.5700,7
13,removelargest,200.0,100.0,20.0,0.9249,14767,0.9901,0.0629,12631.3667,0.0000e+00,0.7982,146.5700,7
13,all,100.0,100.0,20.0,0.7337,8813,0.5755,0.1253,4743.6293,0.0000e+00,0.7050,146.5700,7
13,largestonly,100.0,100.0,20.0,0.8620,5508,0.6570,0.0867,4092.6438,0.0000e+00,0.7050,146.5700,7
13,removelargest,100.0,100.0,20.0,0.8746,7994,0.5800,0.0824,6115.1257,0.0000e+00,0.7050,146.5700,7

17,all,1000.0,100.0,20.0,0.1727,13723,0.3578,0.2983,409.5116,1.4171e-178,0.2925,278.3200,7
17,largestonly,1000.0,100.0,20.0,0.2074,9740,0.4879,0.2823,419.1494,9.2416e-183,0.2925,278.3200,7
17,removelargest,1000.0,100.0,20.0,0.1775,10830,0.4205,0.2959,341.2249,6.4256e-149,0.2925,278.3200,7
17,all,750.0,100.0,20.0,0.4003,14051,0.6936,0.2154,2251.3935,0.0000e+00,0.5285,278.3200,7
17,largestonly,750.0,100.0,20.0,0.4389,9737,0.7977,0.2042,1876.0223,0.0000e+00,0.5285,278.3200,7
17,removelargest,750.0,100.0,20.0,0.4121,10971,0.7427,0.2119,1862.8166,0.0000e+00,0.5285,278.3200,7
17,all,500.0,100.0,20.0,0.6781,13466,0.9790,0.1403,6192.6715,0.0000e+00,0.7067,278.3200,7
17,largestonly,500.0,100.0,20.0,0.6424,10160,0.0665,0.1497,4193.0038,0.0000e+00,0.7067,278.3200,7
17,removelargest,500.0,100.0,20.0,0.6583,11510,0.0258,0.1455,4987.7603,0.0000e+00,0.7067,278.3200,7
17,all,400.0,100.0,20.0,0.7526,15045,0.6853,0.1200,8521.5208,0.0000e+00,0.7503,278.3200,7
17,largestonly,400.0,100.0,20.0,0.7110,10081,0.7430,0.1315,5095.5911,0.0000e+00,0.7503,278.3200,7
17,removelargest,400.0,100.0,20.0,0.7261,11865,0.7149,0.1273,6256.0767,0.0000e+00,0.7503,278.3200,7
17,all,300.0,100.0,20.0,0.8465,19004,0.3410,0.0919,13618.0410,0.0000e+00,0.7864,278.3200,7
17,largestonly,300.0,100.0,20.0,0.7294,10156,0.4121,0.1264,5403.9568,0.0000e+00,0.7864,278.3200,7
17,removelargest,300.0,100.0,20.0,0.8126,14251,0.3826,0.1025,9409.1893,0.0000e+00,0.7864,278.3200,7
17,all,200.0,100.0,20.0,0.9242,14908,0.9429,0.0632,12734.7773,0.0000e+00,0.7982,278.3200,7
17,largestonly,200.0,100.0,20.0,0.8426,11532,0.0351,0.0932,8186.8024,0.0000e+00,0.7982,278.3200,7
17,removelargest,200.0,100.0,20.0,0.8765,13954,0.9925,0.0817,10719.8722,0.0000e+00,0.7982,278.3200,7
17,all,100.0,100.0,20.0,0.6392,9597,0.5892,0.1506,3921.0614,0.0000e+00,0.7050,278.3200,7
17,largestonly,100.0,100.0,20.0,0.7601,7874,0.6258,0.1179,4549.0788,0.0000e+00,0.7050,278.3200,7
17,removelargest,100.0,100.0,20.0,0.7464,8627,0.5987,0.1217,4805.8601,0.0000e+00,0.7050,278.3200,7

18,all,1000.0,100.0,20.0,0.1618,11398,0.3656,0.3038,298.4069,2.5324e-130,0.2913,273.2200,8
18,largestonly,1000.0,100.0,20.0,0.1916,9115,0.5833,0.2893,334.7856,4.0220e-146,0.2913,273.2200,8
18,removelargest,1000.0,100.0,20.0,0.1725,6151,0.5183,0.2984,183.0857,3.0681e-80,0.2913,273.2200,8
18,all,750.0,100.0,20.0,0.4239,11809,0.7146,0.2085,2121.7055,0.0000e+00,0.5280,273.2200,8
18,largestonly,750.0,100.0,20.0,0.4181,9121,0.8620,0.2102,1594.3564,0.0000e+00,0.5280,273.2200,8
18,removelargest,750.0,100.0,20.0,0.4204,6491,0.8268,0.2095,1146.9836,0.0000e+00,0.5280,273.2200,8
18,all,500.0,100.0,20.0,0.7057,11995,0.0025,0.1329,5973.8902,0.0000e+00,0.7066,273.2200,8
18,largestonly,500.0,100.0,20.0,0.6134,9290,0.1059,0.1574,3495.0813,0.0000e+00,0.7066,273.2200,8
18,removelargest,500.0,100.0,20.0,0.7282,8278,0.0866,0.1267,4390.1854,0.0000e+00,0.7066,273.2200,8
18,all,400.0,100.0,20.0,0.7963,13732,0.6957,0.1074,8707.3501,0.0000e+00,0.7500,273.2200,8
18,largestonly,400.0,100.0,20.0,0.7058,9442,0.7770,0.1329,4703.5281,0.0000e+00,0.7500,273.2200,8
18,removelargest,400.0,100.0,20.0,0.7972,9596,0.7577,0.1072,6098.4919,0.0000e+00,0.7500,273.2200,8
18,all,300.0,100.0,20.0,0.8508,18552,0.3624,0.0905,13428.7285,0.0000e+00,0.7870,273.2200,8
18,largestonly,300.0,100.0,20.0,0.7430,9183,0.4319,0.1227,5069.1089,0.0000e+00,0.7870,273.2200,8
18,removelargest,300.0,100.0,20.0,0.8410,12503,0.4191,0.0937,8842.7541,0.0000e+00,0.7870,273.2200,8
18,all,200.0,100.0,20.0,0.9374,14939,0.9557,0.0572,13126.3408,0.0000e+00,0.7986,273.2200,8
18,largestonly,200.0,100.0,20.0,0.8287,10760,0.0644,0.0976,7389.6495,0.0000e+00,0.7986,273.2200,8
18,removelargest,200.0,100.0,20.0,0.9122,14319,0.0136,0.0682,11914.2387,0.0000e+00,0.7986,273.2200,8
18,all,100.0,100.0,20.0,0.7563,8750,0.5794,0.1189,5005.5081,0.0000e+00,0.7050,273.2200,8
18,largestonly,100.0,100.0,20.0,0.8124,7609,0.6289,0.1026,5021.5925,0.0000e+00,0.7050,273.2200,8
18,removelargest,100.0,100.0,20.0,0.9332,7658,0.5880,0.0592,6669.1723,0.0000e+00,0.7050,273.2200,8

30,all,1000.0,100.0,20.0,0.1138,10290,0.6208,0.3318,133.2150,1.3978e-58,0.2913,172.1400,9
30,largestonly,1000.0,100.0,20.0,0.1083,2504,0.0146,0.3356,29.3761,1.7464e-13,0.2913,172.1400,9
30,removelargest,1000.0,100.0,20.0,0.0958,5889,0.7454,0.3447,54.0769,3.2711e-24,0.2913,172.1400,9
30,all,750.0,100.0,20.0,0.3323,10548,0.8977,0.2362,1164.8137,0.0000e+00,0.5275,172.1400,9
30,largestonly,750.0,100.0,20.0,0.2646,2627,0.1553,0.2595,183.9689,1.2685e-80,0.5275,172.1400,9
30,removelargest,750.0,100.0,20.0,0.3482,6209,0.9893,0.2312,752.9756,0.0000e+00,0.5275,172.1400,9
30,all,500.0,100.0,20.0,0.6379,12434,0.1344,0.1509,5059.8034,0.0000e+00,0.7063,172.1400,9
30,largestonly,500.0,100.0,20.0,0.5301,2295,0.3156,0.1793,644.9216,8.2056e-281,0.7063,172.1400,9
30,removelargest,500.0,100.0,20.0,0.6576,8429,0.1933,0.1457,3644.9556,0.0000e+00,0.7063,172.1400,9
30,all,400.0,100.0,20.0,0.7405,13088,0.7903,0.1234,7176.6248,0.0000e+00,0.7501,172.1400,9
30,largestonly,400.0,100.0,20.0,0.6491,2104,0.9492,0.1480,886.5489,0.0000e+00,0.7501,172.1400,9
30,removelargest,400.0,100.0,20.0,0.7510,10025,0.8430,0.1205,5653.3777,0.0000e+00,0.7501,172.1400,9
30,all,300.0,100.0,20.0,0.8216,16718,0.4397,0.0998,11284.3827,0.0000e+00,0.7869,172.1400,9
30,largestonly,300.0,100.0,20.0,0.7394,2095,0.5653,0.1237,1145.3033,0.0000e+00,0.7869,172.1400,9
30,removelargest,300.0,100.0,20.0,0.8045,11968,0.4754,0.1050,7746.0542,0.0000e+00,0.7869,172.1400,9
30,all,200.0,100.0,20.0,0.9487,14966,0.9852,0.0516,13470.6300,0.0000e+00,0.7988,172.1400,9
30,largestonly,200.0,100.0,20.0,0.7778,2568,0.1662,0.1128,1553.4501,0.0000e+00,0.7988,172.1400,9
30,removelargest,200.0,100.0,20.0,0.9242,14629,0.0366,0.0632,12495.1789,0.0000e+00,0.7988,172.1400,9
30,all,100.0,100.0,20.0,0.8313,8230,0.5813,0.0967,5687.6253,0.0000e+00,0.7054,172.1400,9
30,largestonly,100.0,100.0,20.0,0.8731,3818,0.6850,0.0829,2910.7500,0.0000e+00,0.7054,172.1400,9
30,removelargest,100.0,100.0,20.0,0.9607,7597,0.5927,0.0451,7011.1321,0.0000e+00,0.7054,172.1400,9
"""

if __name__ == '__main__':
    pass

