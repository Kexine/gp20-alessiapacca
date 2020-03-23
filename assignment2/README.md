# Assignment 2

Edit this 'README.md' file to report all your results. There is no need to write lengthy reports, just show the requested outputs and screenshots and quickly summarize your observations. Please add your additional files or notes in the folder 'assignment2/results' and refer to or directly show them in this page.

## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.

![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/cat.png)

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).
This is the grid of luigi.off, using PCA in order to align it. <br><br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/luigi.png)
<br>
This is instead the grid of cat.off:
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/catgrid.png)

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results).
The best results were obtained with the simplest shape, the sphere: <br>
By trying with polyDegree = 0, we get <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly0_0.png) <br>
By trying with polyDegree = 1, we get <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly1_1.png) <br>
By trying with polyDegree = 2, we get <br>
![alt text](https://github.com/eth-igl/gp20-alessiapacca/blob/master/assignment2/results/sphere_poly2_2.png) <br>


### Theory question: Save your notes to assignment2/results and add a link to this page.

1) Save your notes and add a link to this page.

2) Show screenshots comparing the 'hound.off' of the normal based reconstruction to the point based reconstruction of the mandatory task.

3) Compare your MLS reconstruction results to the surfaces obtained with Screened Poisson Reconstruction and RIMLS, and try to understand the differences. Report your findings.
