set offset 1,1,1,1
flabel(y)=sprintf("%.2f", y)
plot 'output.csv' using 1:2:(flabel($3)) with labels point offset character 0,character 1 tc rgb "blue"

pause -1