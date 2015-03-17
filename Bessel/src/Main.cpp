#include <stdio.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>
#include <sstream>
#include <fstream>

#define BIN true

using namespace std;


// should be called ./calc_bessel lmin lmax xmin xmax xstep
int main(int argc, char *argv[])
{

    int lmin, lmax;
    float xmin, xmax, x, xstep;
    if (argc == 6){
        istringstream s1(argv[1]), s2(argv[2]), s3(argv[3]), s4(argv[4]), s5(argv[5]);
        s1 >> lmin;
        s2 >> lmax;
        s3 >> xmin;
        s4 >> xmax;
        s5 >> xstep;
    } else {
        cout << "Assuming default values. Run as >> ./calc_bessel lmin lmax xmin xmax xstep << to specify parameters." << endl;
        lmin = 0;
        lmax = 100;
        xmin = 1.0;
        xmax = 200000.0;
        xstep = 1.0;
    }

    if (BIN) {

        ofstream output_binary;
        float bess;
        output_binary.open("bessel_table.bin", ios_base::binary);
        output_binary.write( (char*)&lmin , sizeof(int) );
        output_binary.write( (char*)&lmax , sizeof(int) );
        output_binary.write( (char*)&xmin , sizeof(float) );
        output_binary.write( (char*)&xmax , sizeof(float) );
        output_binary.write( (char*)&xstep , sizeof(float) );
        for (int l = lmin; l <= lmax; ++l) {
            x = xmin;
            cout << "Now writing bessel function for index l = " << l << " ..." << endl;
            while (x <= xmax) {
                bess = boost::math::sph_bessel(l,x);
                output_binary.write((char*)&bess,sizeof(float));
                x += xstep;
            }
        }
        output_binary.close();

    }
    else
    {

        ofstream output;
        output.open("bessel_table.dat");
        output << lmin << " " << lmax << " " << xmin << " " << xmax << " " << xstep << endl; 
        for (int l = lmin; l <= lmax; ++l) {
            x = xmin;
            cout << "Now writing bessel function for index l = " << l << " ..." << endl;
            while (x <= xmax) {
                output << boost::math::sph_bessel(l,x) << " ";
                x += xstep;
            }
            output << endl;
        }
        output.close();
    }

    return 0; 
}
