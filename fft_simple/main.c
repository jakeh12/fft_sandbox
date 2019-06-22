//
//  main.c
//  fft_simple
//
//  Created by Jakub Hladik on 6/22/19.
//  Copyright © 2019 Jakub Hladik. All rights reserved.
//

#include <stdio.h>
#include <math.h>


void twiddle_factor_lookup_table(int n, double twiddle_lut_r[], double twiddle_lut_i[])
{
    double theta;
    int i;
    for (i = 0; i < n/2; i++) {
        theta = ((-2 * M_PI) / n) * i;
        twiddle_lut_r[i] = cos(theta);
        twiddle_lut_i[i] = sin(theta);
    }
}

void butterfly(double x0_r, double x0_i, double x1_r, double x1_i, double w_r, double w_i) {
    //printf("inputs:\nx0: %.4f%+.4fi\tx1: %.4f%+.4fi\tw: %.4f%+.4fi\n", x0_r, x0_i, x1_r, x1_i, w_r, w_i);
    double y0_r, y0_i, y1_r, y1_i, x0_r_minus_x1_r, x0_i_minus_x1_i;
    // calculate y0 (complex y0 = x0 + x1)
    y0_r = x0_r + x1_r;
    y0_i = x0_i + x1_i;
    // calculate y1 (complex y1 = w * (x0 - x1))
    x0_r_minus_x1_r = x0_r - x1_r;
    x0_i_minus_x1_i = x0_i - x1_i;
    y1_r = w_r * x0_r_minus_x1_r - w_i * x0_i_minus_x1_i;
    y1_i = w_r * x0_i_minus_x1_i + w_i * x0_r_minus_x1_r;
    //printf("outputs:\ny0: %.4f%+.4fi\ty1: %.4f%+.4fi\n\n", y0_r, y0_i, y1_r, y1_i);
}

int main(int argc, const char * argv[]) {
    int n = 8; // n has to be a power of two
    double twiddle_lut_r[n], twiddle_lut_i[n];
    twiddle_factor_lookup_table(n, twiddle_lut_r, twiddle_lut_i);
    int i;
    for (i = 0; i < n/2; i++) {
            printf ("%+.10f%+.10fi\n", twiddle_lut_r[i], twiddle_lut_i[i]);
    }
    
    butterfly(1.224242, 0.113198, 5.422113, -0.411351, 0.00113141, -0.948824);
    return 0;
}
