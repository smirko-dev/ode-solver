set terminal png
set datafile separator ','
set output "lorenz.png"
plot "lorenz.csv" using 1:2 with lines linecolor 'blue'
unset output

