# ShadowLineDistributions


Code and data for the paper ``Shadow line distributions'' (Jennifer S. Balakrishnan, Mirela Çiperiani, Barry Mazur, Karl Rubin), 
https://arxiv.org/abs/2409.00891

Each shadow line distribution is computed using the corresponding .sage file in its subfolder (which makes a call to shadow.sage, currently in the top level of the directory) . 
The output data is recorded in its subfolder. 

For example, to produce the shadow line distribution for (433.a1, 3), the file 433.a1_3.sage in the subfolder anomalous/433.a1_3 sets up the computation, and output should match the data recorded in anomalous/433.a1_3/433.a1_3_out.


To compute 3-adic heights, replace the following Sage library file:
* src/sage/schemes/elliptic_curves/padics.py
  
with the file padics.py in this repository and rebuild Sage. This is built on Sage-9.8.
