set term pdfcairo
set size ratio -1
set key off
set xrange [-5.55112e-17:4.7124]
set yrange [-1:1]
set style line 1 lw 1.5 lt 1 pt 7 ps 0.5
$points << EOD
0.3927 0.5
1.1781 -0.5
1.9635 0.5
2.7489 -0.5
3.5343 0.5
4.3197 -0.5
EOD
plot sample [0.3927:1.1781] (0 * sinh(10 * (1.1781 - x)) + 15.8202 * sinh(10 * (x - 0.3927))) / 128801 + 0.5 * (1.1781 - x) / 0.7854 + -0.658202 * (x - 0.3927) / 0.7854 ls 1, \
[1.1781:1.9635] (15.8202 * sinh(10 * (1.9635 - x)) + -16.9673 * sinh(10 * (x - 1.1781))) / 128801 + -0.658202 * (1.9635 - x) / 0.7854 + 0.669673 * (x - 1.1781) / 0.7854 ls 1, \
[1.9635:2.7489] (-16.9673 * sinh(10 * (2.7489 - x)) + 16.9673 * sinh(10 * (x - 1.9635))) / 128801 + 0.669673 * (2.7489 - x) / 0.7854 + -0.669673 * (x - 1.9635) / 0.7854 ls 1, \
[2.7489:3.5343] (16.9673 * sinh(10 * (3.5343 - x)) + -15.8202 * sinh(10 * (x - 2.7489))) / 128801 + -0.669673 * (3.5343 - x) / 0.7854 + 0.658202 * (x - 2.7489) / 0.7854 ls 1, \
[3.5343:4.3197] (-15.8202 * sinh(10 * (4.3197 - x)) + 0 * sinh(10 * (x - 3.5343))) / 128801 + 0.658202 * (4.3197 - x) / 0.7854 + -0.5 * (x - 3.5343) / 0.7854 ls 1, \
$points with points ls 1