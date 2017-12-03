# Extended Kalman Filter
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

This project uses an extended kalman filter to continuously estimate the state of a moving vehicle in a simulator. The vehicle reports noisy lidar and radar data continuously.

The simulator used can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

## Building the project

To build the project simply run the build.sh script in the root directory of the project, ensure that all the listed dependencies are installed on your system

`
$ ./build.sh
`

## Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Code Style

Stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).
