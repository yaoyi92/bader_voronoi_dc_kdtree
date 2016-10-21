# bader_voronoi_dc_kdtree

voronoi charge analysis is sometimes used in my research. The bader code is not optimized for this functionality (http://theory.cm.utexas.edu/henkelman/code/bader/).
I borrowed some algorithms from computational science to accelarate the voronoi charge separation. 
It is now 100s times faster and only takes several seconds for a cube file used to cost half an hour for the analysis.

for voronoi analysis

example:
bader -n bader -c voronoi test.cube

===
speed test (32 h2o system)
===
time spend (s) for voronoi analysis with different method.
mesh size  original code dc code     dc_kdtree code
100^3       42.32          4.83         0.37
144^3       126.31         11.14        0.85
192^3       299.97         19.73        1.52
270^3       833.57         44.67        3.42

===
speed test (512 Si + 1 H)
mesh size  original code dc code     dc_kdtree code
192^3       1713.11        184.98       2.23

With dc_kdtree code, loading the cube file takes more time ;)

===
method
===
1. divide and conquer for the mesh of electrons
2. kdtree for fast closest atom search ( expand the atoms to periodic images to take account of PBC)
   https://github.com/jmhodges/kdtree2

YY
10.2016
