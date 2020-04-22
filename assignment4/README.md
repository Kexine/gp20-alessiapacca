# Assignment 4

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment4/results' and refer to or directly show them in this page.


## Required results

* Screenshots of the parameterizations and textured (checkerboard) models for all the implemented methods and boundary conditions (models: cathead.obj, hemisphere.off, hemisphere_non_convex_boundary.off,Octo_cut2.obj)

<br>
Let's see the parametrizations with the shape 'cat.off'.<br>
By displaying the uniform laplacian, we get:

![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat1.png) <br><br>
By displaying the cotangent laplacian, which is the same of LSCM in case of constrained boundary, we get:<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat2_3.png) <br><br>
I also computed the cotangent Laplacian without using the function 'igl::cotmatrix', and the result did not changed: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat2_3_cotlap.png) <br><br>
The parametrization with ARAP with constrained boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4.png) <br><br>
The parametrization with ARAP with free boundary outputs a better results: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4_freeboundary%20.png) <br><br>
<br>





* Several examples of the distortion visualizations.
Let's see the distorsions with the same model showed before. <br>
Let's start from the shape cat.off'.<br>
We can visualize the distorsion with angle preservation for the uniform laplacian: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat1_distA.png) <br><br>
Again, we can visualize the distorsion with angle preservation for cotangent laplacian and for LSCM, in case of constrained boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat2_3_distA.png) <br><br>
The distorsion with length preservation instead, with ARAP parametrization is: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4_freebound_distL.png)


