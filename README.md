# Extended Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

---

## Results

### Dataset 1

![data/output-1.png](Output for dataset 1)

### Dataset 2

![data/output-2.png](Output for dataset 1)

## Notes

I made a few changes after getting the original code working, as detailed below.

### Avoid explicitly inverting S to find the Kalman gain

It also seems to improve performance by 5%-10% on dataset 1, even though the
matrices are very small.

### Use fixed-size matrix types

Eigen provides both fixed and dynamic-sized matrices; fixed size matrices
provide more compile-time correctness checking and should in principle be
faster by allowing the compiler to calculate more at compile time.

Refactoring to use fixed sized matrices turned out to require some fairly large
changes, because it prevented me from using the technique of switching out the H
and R matrices in the filter depending on which sensor measurement we're using.

Instead, it required extracting the update step into two separate (templated)
Sensor classes, each with H and R matrices of the appropriate size. It also
required moving the Jacobian calculation out of Tools and into the Radar Sensor
class.

The performance benefits of doing this turned out to be pretty small: about
another 5% on dataset 1.

However, architecturally, it seems much clearer, and (subjectively) it cleaned
up the code a lot.

### Make better use of first measurement

In the lecture, we used the position components of the first measurement to initialise the state and hardcoded the initial position uncertainties at 1 and the initial velocity uncertainties at 1000. This works fine on the sample data, but it seems risky to take the first measurement without any filtering. I tried two alternative approaches.

### Set a prior and treat the first measurement as a normal update

Given that we know that the obstacle is a pedestrian, and that our sensors have limited range, we can set a reasonable prior on the position and speed. For example, if we set the position variance to 1000 and the velocity variance to 25, our prior is that the pedestrian is within ~50m and moving at <8m/s in either x or y, which seems reasonable.

However, this did not work, for several reasons:

   - Setting the prior position mean to zero, which is the only choice that doesn't add bias in a particular direction, causes us to skip the first update if it is a radar update, because the `h` function is undefined when the range is zero. It was therefore necessary to set the position to a small value, such as (0.1, 0.1), which gives a preference for the first quadrant.

   - The Taylor expansions in the linearised `Hj` matrix are expanded about the position from the current state. On a normal update, the state has been updated by the prediction step, but before the first measurement, the state values are zero (or close to zero) from the prior. If the first measurement of position is far from the origin, then the expansion for `Hj` can be very inaccurate.

To work around both problems, we can take the position from the first measurement, use it to initialise `Hj`, and then run an otherwise-normal update. This approach seems to work quite well. For example, the first measurement in dataset 1 is
```
rho=8.46642 phi=0.0287602 rho_dot=-3.04035
```
which means that the obstacle is near the positive x axis and moving toward the origin at 3m/s. If we use `rho` and `phi` to find the position (and set the initial `vx` and `vy` to zero) to use in the calculation of the initial `Hj`, the state covariance matrix (`P`) after updating with the first measurement is:
```
  0.0899708  0.000732516           0           0
0.000732516    0.0645292           0           0
          0            0    0.110276   -0.716031
          0            0   -0.716031     24.9794
```
The lower right block gives a Gaussian with little uncertainty in the x direction (constrained by the measurement) and essentially unchanged uncertainty in the y direction (the prior velocity variance was 25). The negative covariance terms tilt the Gaussian slightly, perpendicular to the radial with angle `phi`.

However, this makes the code more complicated, the effects are only significant in the first 3-4 iterations, and it may be deviating too far from the spec for this project. I've therefore left it alone.

### Set initial position uncertainties based on measurement uncertainties

For the initial state position uncertainties, it seems more principled to use the measurement uncertainty matrix to set the initial process noise, rather than using arbitrary values.

For the laser, which measures the `px` and `py` components of the state directly, we can just use the relevant entries from the `R` matrix for the laser. A single laser measurement provides no data on speed, so we just leave the initial velocity variances at an arbitrary value of 1000.

For the radar, we can approximate the uncertainty in the inferred `px` and `py` positions based on the uncertainties in the `rho` and `phi` measurements by using the usual approximations for [propagation of uncertainty](https://en.wikipedia.org/wiki/Propagation_of_uncertainty) through a product and the cosine and sine functions.

The radar does provide some information about speed, but it is only the speed in the radial direction of the object (from Doppler). To decompose the radial speed into `vx` and `vy` components, we'd have to invert the `h` function, but unfortunately it is not invertible. The simplest way that I have found to incorporate the speed information from the radar is the approach in the section above --- while `h` is not invertible, the linearised mapping defined by `Hj` is, after transformation in the filter, invertible.

Compared to the original RMSEs with hard coded uncertainties, these changes to the initialisation have very small and mixed effects on the final RMSEs.

The final RMSEs for dataset 1 (which has a radar measurement as its first measurement):

```
final RMSE(px): 0.0651649 -> 0.0651493
final RMSE(py): 0.0605378 -> 0.0605818
final RMSE(vx): 0.54319   -> 0.543658
final RMSE(vy): 0.544191  -> 0.54431
```

The final RMSEs for dataset 2 (which has a laser measurement as its first measurement):

```
final RMSE(px): 0.185496 -> 0.185514
final RMSE(py): 0.190302 -> 0.190295
final RMSE(vx): 0.476754 -> 0.475961
final RMSE(vy): 0.804469 -> 0.805237
```

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./ExtendedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./ExtendedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!

More information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/12dd29d8-2755-4b1b-8e03-e8f16796bea8)
for instructions and the project rubric.
