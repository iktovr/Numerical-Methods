set term pdfcairo
set size ratio -1
set key off
set xrange [0.717:9.753]
set yrange [-0.5:5.5]
set style line 1 lw 1.5 lt 1 pt 7 ps 0.5
$points << EOD
1.47 5
2.22 3.68
2.6 1.38
3.35 0.53
4.5 0.1
4.86 0.17
5.52 2
6 2.96
7.5 3.68
9 3.92
EOD
plot sample [1.47:2.22] (0 * sinh(0.02 * (2.22 - x)) + -15.1316 * sinh(0.02 * (x - 1.47))) / 6.00023e-06 + 5 * (2.22 - x) / 0.75 + 37832.7 * (x - 1.47) / 0.75 ls 1, \
[2.22:2.6] (-15.1316 * sinh(0.02 * (2.6 - x)) + 22.2139 * sinh(0.02 * (x - 2.22))) / 3.04003e-06 + 37832.7 * (2.6 - x) / 0.38 + -55533.4 * (x - 2.22) / 0.38 ls 1, \
[2.6:3.35] (22.2139 * sinh(5.46174 * (3.35 - x)) + -0.328081 * sinh(5.46174 * (x - 2.6))) / 896.431 + 0.635333 * (3.35 - x) / 0.75 + 0.540998 * (x - 2.6) / 0.75 ls 1, \
[3.35:4.5] (-0.328081 * sinh(4.41076 * (4.5 - x)) + 0.168136 * sinh(4.41076 * (x - 3.35))) / 1551.96 + 0.546864 * (4.5 - x) / 1.15 + 0.0913576 * (x - 3.35) / 1.15 ls 1, \
[4.5:4.86] (0.168136 * sinh(0.02 * (4.86 - x)) + 8.85515 * sinh(0.02 * (x - 4.5))) / 2.88002e-06 + -420.239 * (4.86 - x) / 0.36 + -22137.7 * (x - 4.5) / 0.36 ls 1, \
[4.86:5.52] (8.85515 * sinh(0.02 * (5.52 - x)) + -4.02308 * sinh(0.02 * (x - 4.86))) / 5.28015e-06 + -22137.7 * (5.52 - x) / 0.66 + 10059.7 * (x - 4.86) / 0.66 ls 1, \
[5.52:6] (-4.02308 * sinh(0.02 * (6 - x)) + -2.72526 * sinh(0.02 * (x - 5.52))) / 3.84006e-06 + 10059.7 * (6 - x) / 0.48 + 6816.11 * (x - 5.52) / 0.48 ls 1, \
[6:7.5] (-2.72526 * sinh(2.77653 * (7.5 - x)) + -0.209643 * sinh(2.77653 * (x - 6))) / 248.093 + 3.31351 * (7.5 - x) / 1.5 + 3.70719 * (x - 6) / 1.5 ls 1, \
[7.5:9] (-0.209643 * sinh(2.77653 * (9 - x)) + 0 * sinh(2.77653 * (x - 7.5))) / 248.093 + 3.70719 * (9 - x) / 1.5 + 3.92 * (x - 7.5) / 1.5 ls 1, \
$points with points ls 1