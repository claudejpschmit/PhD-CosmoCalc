PhD-CosmoCalc
=============

Contains Cosmology calculator ala Ned Wright.

Setup
-----

Download and unzip into a working directory. Then cd into working directory.
The setup for this program requires the setup for 3 different dependencies.

- Setup for CAMB: 

```sh
$ cd CAMB/
```

view readme.html for CAMB if needed at this point. Compile CAMB and test run it:

```sh
$ make
$ ./camb params.ini
```

The Makefile for CAMB may have to be eddited depending on the Fortran compiler present on your machine.
The test call of camb is important as it writes a binary file of spherical bessels read in by the CosmoCalculator.

- Setup for Bessel function wrapper

cd back to working directory, then:

```sh
$ cd Bessel/Bessel_from_CAMB/
$ f2py -c -m bessels test_bessel.F90
```

This creates the object bessels.so which is called by the CosmoCalculator.

- Setup for camb4py

cd back to working directory, then:

```sh
$ cd camb4py/
$ python setup.py build --no-builtin
```

This generates a build/ folder which is used by the CosmoCalculator.

Usage
-----

If the setup has been correct you should be able to cd back to the working directory and call

```sh
$ python calc.py 
$ python plot.py 
$ python writer.py 
```

add -h to any of the above to gain more information 
