set title 'Poission equation solution - Dirichlet boundary conditions-Explicit'
set key off
#set palette grey
#set pm3d map
set pm3d
set ylabel 'y'
set xlabel 'x'
set zlabel 'T(x,y)'
set hidden3d
splot 'Grid.dat' matrix
replot



