# caustics_via_refraction
Gradient-based optimisation method for discrete refractive surface generation, for goal-based caustics.
Modified implementation of the paper "Poisson-Based Continuous Surface Generation for Goal-Based Caustics", by Yue et al.(http://www.cs.columbia.edu/~yonghao/tog14/absttog14.html). In this version I use gradient descent optimization method, instead of the Poisson based, used in the paper.
I also don't perform a second optimization in order to obtain a continuos surface.
The results are quite satisfing with simple gradient input images, as the following:

![test3](https://user-images.githubusercontent.com/77103965/178326794-80bc4d1c-f7f1-4264-aff3-e6184dd5ead4.png)
![visualization30](https://user-images.githubusercontent.com/77103965/178326894-995821dd-2210-4271-9763-8432437dd703.png)

The visualization is based on corresponding deformed area ratios, but the normal are also computed in order to implement a physical-based simulation on blender.
