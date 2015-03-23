# SETUP script to build dependencies of the CosmoCalculator
# Building CAMB
cd CAMB/
make

# Running basic test for CAMB
./camb params.ini

# Installing F2Py
cd ../F2Py/
python setup.py install

# Building Spherical Bessel function
cd ../Bessel/Bessels_from_CAMB/
f2py -c -m bessels test_bessel.F90

# Building camb4py wrapper
cd ../../camb4py-master/
python setup.py build --no-builtin

echo "The Setup was successful if no errors arose. Test by typing: python plot.py"
