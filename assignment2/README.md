# Assignment 2

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment2/results' and refer to or directly show them in this page.

## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat.png)

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).
This is the grid of luigi.off, using PCA in order to align it. <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/luigi.png) <br><br>
This is instead the grid of cat.off: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/catgrid.png) <br>

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results). <br><br>
The best results were obtained with the simplest shape, the sphere. <br>
I used resolution = 30 for all of the axis, wendLandRadius rate = 0.1, and initial epsilon = 0.03. WendLandRadius rate was multiplied for the diagonal of the bounding box, in order to be adapted to every shape. <br>
By trying with polyDegree = 0, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly0_0.png) <br><br>
By trying with polyDegree = 1, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly1_1.png) <br><br>
By trying with polyDegree = 2, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly2_2.png) <br><br><br>
We can see that with higher degree, we have some artifacts, that are more visible in other more complex shapes. <br>
We can see the example of the bunny. <br>
By trying with polyDegree = 0, we get: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_0.png) <br><br>
I also tried to use anisotropic resolutions for the axis, but it did not change the result in a relevant way: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_1.png) <br><br>
By trying with polyDegree = 1, we already get some artifacts, compared to the sphere's result: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_2.png) <br><br>
By trying with polyDegree = 2, I had to change a bit the values of the wendLandRadius and of the initial epsilon: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_3.png) <br><br>
Here are the results by using bunny_1000,horse and hound.off with polyDegree = 0: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/bunny_500_0.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/horse_0.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/hound_0.png) <br><br><br>
I also experimented with the cat shape.<br>
First I tried polyDegree = 0 with different resolutions of the cat: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly0_0.png) <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly0_1.png) <br>
With higher resolutions, we are able to have smoother shapes. <br>
By using polyDegree = 1, the result is: <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly1_2.png) <br><br>
The artifacts can be limited if we increase the wendLandRadius rate, obtaining a "blob" shape: 
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat_poly1_3.png) <br><br>



Theory question


### Theory question: Save your notes to assignment2/results and add a link to this page.

1) Save your notes and add a link to this page.

2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.

3) Compare your MLS reconstruction results to the surfaces obtained with Screened Poisson Reconstruction and RIMLS, and try to understand the differences. Report your findings.

4) I ---descrizione metodo--. 
We can see differences in Luigi shape: 
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/collage.jpg) <br><br>

