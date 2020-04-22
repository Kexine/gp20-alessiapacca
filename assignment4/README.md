# Assignment 4

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment4/results' and refer to or directly show them in this page.


## Required results

* Screenshots of the parameterizations and textured (checkerboard) models for all the implemented methods and boundary conditions (models: cathead.obj, hemisphere.off, hemisphere_non_convex_boundary.off,Octo_cut2.obj)

All the screenshots, for all the models, can be found in the folder [results](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results).
<br>
Here I am gonna show some results for two of the models, the `cat.off` and the `hemisphere`. 

<br>
Let's see the parametrizations with the shape `cat.off`.<br>
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

Let's now show the parametrizations with the `hemisphere` shape: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem1.png) <br><br>
By displaying the cotangent laplacian, which is the same of LSCM in case of constrained boundary, we get:<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem2_3.png) <br><br>
Instead, in the case of free boundary, LSCM is different: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem3_freebound.png) <br><br>
The parametrization with ARAP with free boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem4_freebound.png) <br><br>



* Several examples of the distortion visualizations.
Let's see the distorsions with the same model showed before. <br>
Let's start from the shape `cat.off`.<br>
We can visualize the distorsion with angle preservation for the uniform laplacian: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat1_distA.png) <br><br>
Again, we can visualize the distorsion with angle preservation for cotangent laplacian and for LSCM, in case of constrained boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat2_3_distA.png) <br><br>
The distorsion with length preservation instead, with ARAP parametrization is: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4_freebound_distL.png)

Let's now see the `hemisphere`.<br>
We can visualize the distorsion with angle preservation for the uniform laplacian: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_1_distA.png) <br><br>
Again, we can visualize the distorsion with angle preservation for cotangent laplacian and for LSCM, in case of constrained boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_2_distA.png) <br><br>
The distorsion with angle preservation with ARAP parametrization is: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_4_freebound_distA.png)


