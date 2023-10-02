# approximations
Implementations of sin, cos, atan, exp, and log in ANSI C, with speed and accuracy similar to glibc.

![Graphs of the different functions.](https://github.com/calebabutler/approximations/blob/main/graphs.png?raw=true)

How to Test
===========

To test the functions, you need git, GCC, GNU make, glibc, GNU time, python3,
matplotlib, and numpy. If you are on a Debian or Ubuntu machine, you can just
run the command:

    > sudo apt-get install build-essential git time python3 python3-matplotlib python3-numpy

The, run the following commands to test:

    > git clone https://github.com/calebabutler/approximations.git
    > cd approximations
    > make test

