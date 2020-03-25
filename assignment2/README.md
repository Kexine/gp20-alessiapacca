# Assignment 2

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment2/results' and refer to or directly show them in this page.

## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat.png)

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).<br>
This is the grid of luigi.off, using PCA in order to align it. <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/luigi.png) <br><br>
This is instead the grid of cat.off: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/catgrid.png) <br>

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results). <br><br>
The best results were obtained with the simplest shape, the sphere. <br>
I used resolution = 30 for all of the axis, wendLandRadius rate = 0.1, and initial epsilon = 0.03. WendLandRadius rate was multiplied for the diagonal of the bounding box, in order to be adapted to every shape. <br>
By trying with ```polyDegree``` = 0, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly0_0.png) <br><br>
By trying with ```polyDegree``` = 1, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly1_1.png) <br><br>
By trying with ```polyDegree``` = 2, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly2_2.png) <br><br><br>
We can see that with higher degree, we have some artifacts, that are more visible in other more complex shapes. <br>
We can see the example of the bunny. <br>
By trying with ```polyDegree``` = 0, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_0.png) <br><br>
I also tried to use anisotropic resolutions for the axis, which gives us a "stretched" shape: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_1.png) <br><br>
By trying with ```polyDegree``` = 1, we already get some artifacts, compared to the sphere's result: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_2.png) <br><br>
By trying with ```polyDegree``` = 2, changing the values of wendLandRadius Rate and of the initial epsilon: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_3.png) <br><br>
Here are the results by using ```bunny_1000```,```horse``` and ```hound.off``` with ```polyDegree``` = 0: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_0.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/horse_0.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/horse_0_1.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/hound_0.png) <br><br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/hound_0_1.jpg) <br><br><br>
I also experimented with the cat shape.<br>
First I tried ```polyDegree``` = 0 with different resolutions of the cat: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly0_0.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly0_1.png) <br>
With higher resolutions, we are able to have smoother shapes. <br>
By using ```polyDegree``` = 1, the result is: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly1_2.png) <br><br>
The artifacts can be limited if we increase the wendLandRadius rate, obtaining a "blob" shape: 
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly1_3.png) <br><br>



### Theory question: Save your notes to assignment2/results and add a link to this page.

1) Save your notes and add a link to this page.<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/es1.jpg) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/es2.jpg) <br><br>


2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.<br>
The method, instead of using MLS to blend between the constrained values linked to each polygon/point, blends between functions associated with them. It does that by updating the vector of the values of the neighboring point cloud points (not the constrained ones), on the right part of the system we solve, by a factor that depends on the normal of the point in consideration. 
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/screen.png) <br>
One advantage of this approach is not having to use offset constraints, which is faster when solving. Moreover, using the additional constraints in MLS cause an oscillatory behaviour when the evaluation point is away from the surface. The reason for this side effect is that the distance between the evaluation point and the point on the surface is larger than the offset distance, leading the inside and outside contraints to cancel out.<br>
Let's see the comparison between the two methods, with the hound shape.<br>
By using the right parameters, we see a small difference in the neck, in the mouth and in the ears. <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/normalconfronto.jpg) <br><br>
We can also see differences in Luigi shape:<br>
On the left, Luigi reconstructed without the normal constrain; on the right, Luigi reconstructed with the normal constrain.
Some features, like the hat or the arms of the shape, are better reconstructed with the constrain. <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/collage.jpg) <br>
This method exhibits little oscillation, and it's also way more precise. 


3) Compare your MLS reconstruction results to the surfaces obtained with Screened Poisson Reconstruction and RIMLS, and try to understand the differences. Report your findings.<br><br>
Here are some results found with Screened Poisson reconstruction:<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/catpoisson.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/horsepoisson.png) <br>
Screened Poisson, even if used with low resolutions, has an accuracy which exceeds the one of our MLS reconstruction. Moreover, it has significantly faster processing times. Therefore, it can handle large models in less time, resulting in a reconstruction which is much more detailed and precise.<br> As the paper explains, this is probably because of hierarchical clustering of the points and of a conforming octree structure, that enable a multigrid algorithm with linear complexity on the number of input points <br>
Moreover, the accuracy of our MLS method is also influenced by the sensitivity to outliers and the smoothing out of
small features.<br><br>
Another solution to MLS was developed with RIMLS, which is a robust Implicit MLS, that aims to preserve sharp edges and to be robust to outliers. In this case, the results better preserve fine details and handle sharp features. The reason is that this method uses robust kernel regression combined with MLS. <br>
Here are some results found with RIMLS:<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/catmarch.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/horsemarch.png) <br><br>
In particular, we can look at the reconstruction of Luigi with this method. On the left, there is the reconstruction with our method, on the right, the reconstruction with RIMLS: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/luigiconfronto.jpg) <br><br>
By changing some parameters:<br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/luigiconfronto1.jpg) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/luigiconfronto2.jpg) <br><br>
In this case it's clear that this approach preserves the shape in a better way, without instabilities, and this can be seen from the legs which are more separated, or from the nose which is less mashed into the face of Luigi, and from the hat's reconstruction which is more stable. Most important, it's evident that the corners are better preserved with RIMLS, while with our method the shape is more smoothed.
<br>


