# Assignment 4

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment4/results' and refer to or directly show them in this page.


## Required results

* Screenshots of the parameterizations and textured (checkerboard) models for all the implemented methods and boundary conditions (models: cathead.obj, hemisphere.off, hemisphere_non_convex_boundary.off,Octo_cut2.obj)

All the screenshots, for all the models, can be found in the folder [results](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results).
<br>

<br>
Let's see the parametrizations with the shape `cat.off`.<br>
By displaying the uniform laplacian, we get:

![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat1.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat1_par.png) <br><br>
By displaying the cotangent laplacian, which is the same of LSCM in case of constrained boundary, we get:<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat2_3.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat2_par.png) <br><br>
I also computed the cotangent Laplacian without using the function 'igl::cotmatrix', and the result did not changed: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat2_3_cotlap.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat_3_cotlap_par.png) <br><br>
The parametrization with ARAP with constrained boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4_par.png) <br><br>
The parametrization with ARAP with free boundary outputs a better results: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4_freeboundary%20.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cat4_free_par.png) <br><br>
<br>

Let's now show the parametrizations with the `hemisphere_non_convex_boundary.off` shape: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem1.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem1_par.png) <br><br>
By displaying the cotangent laplacian, which is the same of LSCM in case of constrained boundary, we get:<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem2_3.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem2_3_par.png) <br><br>
Instead, in the case of free boundary, LSCM is different: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem3_freebound.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem3_free_par.png) <br><br>
The parametrization with ARAP with free boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem4_freebound.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem4_free_par.png) <br><br>

Results for the `hemisphere.off`:<br>
Key 1:
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_1.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_1_par.png) <br><br>
Key 2 and 3:
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_2.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_2_par.png) <br><br>
Key 3, with free boundary: 
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_3_freebound.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_3_freebound_par.png) <br><br>
Key 4: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_4.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_4_par.png) <br><br>
Key 4, with free boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_4_freebound.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/hem_nonconvex_4_freebound_par.png) <br><br>


Results for the `Octo_cut2.obj`:<br>
Key 1:
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_1.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_1_par.png) <br><br>
Key 2 and 3:
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_2_3.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_2_par.png) <br><br>
Key 3, with free boundary: 
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_3_freebound.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_3_freebound_par.png) <br><br>
Key 4: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_4.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_4_par.png) <br><br>
Key 4, with free boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_4_freebound.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/octo_4_freebound_par.png) <br><br>




Results for the `cow`:<br>
Key 1:
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cow_1.png) <br><br>
Key 2 and 3:
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cow_2_3.png) <br><br>
Key 3, with free boundary: 
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cow_3_freebound.png) <br><br>
Key 4: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cow_4.png) <br><br>
Key 4, with free boundary: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment4/results/cow_4_freebound.png) <br><br>


* Several examples of the distortion visualizations.<br>
For the distorsion measure I tried implementing several approaches, as described in the course's slides.<br>
I tried with: <br> 
1. The energies formulas of the two parametrizations. <br>
2. The maximum between the smallest eigenvalue and the reciprocal of the biggest eigenvalue of J. <br>
3. The square root of the sum of the squares of the biggest and the smallest eigenvalues of J. <br><br>
At the end I used the energies as the distorsion measures.<br>
Let's see the results with some of the models showed before. <br>
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


