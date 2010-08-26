plot 'stats.txt' using 1:4 title '0 0 ', \
     'stats.txt' using 1:7 title '1 0 ', \
     'stats.txt' using 1:10 title '0 1 ', \
     'stats.txt' using 1:13 title '1 1 ', \
     'stats.txt' using 1:14 title 'avr error'
set key 25000, 0.5
set terminal postscript eps enhanced color
set output "file.eps"
replot
