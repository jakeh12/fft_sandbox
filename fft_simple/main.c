//
//  main.c
//  fft_simple
//
//  Created by Jakub Hladik on 6/22/19.
//  Copyright © 2019 Jakub Hladik. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <complex.h>

void generate_twiddle_factor_lookup_table(int n, complex double twiddle_lut[])
{
    int r;
    for (r = 0; r < n; r++) {
        twiddle_lut[r] = cexp(((-I * 2 * M_PI) / n) * r);
    }
}

void butterfly(complex double x0, complex double x1, complex double wk, complex double *y0, complex double *y1) {
    *y0 = x0 + wk * x1;
    *y1 = x0 - wk * x1;
}

int main(int argc, const char * argv[]) {
    int n = 8; // n has to be a power of two
    complex double twiddle_lut[n];
    generate_twiddle_factor_lookup_table(n, twiddle_lut);
    printf("twiddle factors:\n");
    int r;
    for (r = 0; r < n; r++) {
        printf("%+f%+fi\n", crealf(twiddle_lut[r]), cimagf(twiddle_lut[r]));
    }
    
    complex double input[n];
    input[0] = +0.0 + 0.0*I;
    input[4] = +0.0 + 0.0*I;
    input[2] = +1.0 + 0.0*I;
    input[6] = +1.0 + 0.0*I;
    input[1] = +0.0 + 0.0*I;
    input[5] = +0.0 + 0.0*I;
    input[3] = +1.0 + 0.0*I;
    input[7] = +1.0 + 0.0*I;
    
    // stage 0
    complex double y0_0_0, y1_0_0;
    butterfly(input[0], input[4], twiddle_lut[0], &y0_0_0, &y1_0_0);//
    complex double y0_0_1, y1_0_1;
    butterfly(input[2], input[6], twiddle_lut[0], &y0_0_1, &y1_0_1);//
    complex double y0_0_2, y1_0_2;
    butterfly(input[1], input[5], twiddle_lut[0], &y0_0_2, &y1_0_2);//
    complex double y0_0_3, y1_0_3;
    butterfly(input[3], input[7], twiddle_lut[0], &y0_0_3, &y1_0_3);//

    // stage 1
    complex double y0_1_0, y1_1_0;
    butterfly(y0_0_0, y0_0_1, twiddle_lut[0], &y0_1_0, &y1_1_0);//
    complex double y0_1_1, y1_1_1;
    butterfly(y1_0_0, y1_0_1, twiddle_lut[2], &y0_1_1, &y1_1_1);//
    complex double y0_1_2, y1_1_2;
    butterfly(y0_0_2, y0_0_3, twiddle_lut[0], &y0_1_2, &y1_1_2);//
    complex double y0_1_3, y1_1_3;
    butterfly(y1_0_2, y1_0_3, twiddle_lut[2], &y0_1_3, &y1_1_3);//
    
    // stage 2
    complex double y0_2_0, y1_2_0;
    butterfly(y0_1_0, y0_1_2, twiddle_lut[0], &y0_2_0, &y1_2_0);//
    complex double y0_2_1, y1_2_1;
    butterfly(y0_1_1, y0_1_3, twiddle_lut[1], &y0_2_1, &y1_2_1);//
    complex double y0_2_2, y1_2_2;
    butterfly(y1_1_0, y1_1_2, twiddle_lut[2], &y0_2_2, &y1_2_2);//
    complex double y0_2_3, y1_2_3;
    butterfly(y1_1_1, y1_1_3, twiddle_lut[3], &y0_2_3, &y1_2_3);//
    
    complex double output[n];
    output[0] = y0_2_0;
    output[1] = y0_2_1;
    output[2] = y0_2_2;
    output[3] = y0_2_3;
    output[4] = y1_2_0;
    output[5] = y1_2_1;
    output[6] = y1_2_2;
    output[7] = y1_2_3;
    
    printf("input vector:\n");
    printf("%+f%+fi\n", crealf(input[0]), cimagf(input[0]));
    printf("%+f%+fi\n", crealf(input[4]), cimagf(input[4]));
    printf("%+f%+fi\n", crealf(input[2]), cimagf(input[2]));
    printf("%+f%+fi\n", crealf(input[6]), cimagf(input[6]));
    printf("%+f%+fi\n", crealf(input[1]), cimagf(input[1]));
    printf("%+f%+fi\n", crealf(input[5]), cimagf(input[5]));
    printf("%+f%+fi\n", crealf(input[3]), cimagf(input[3]));
    printf("%+f%+fi\n", crealf(input[7]), cimagf(input[7]));

    printf("output vector:\n");
    for (r = 0; r < n; r++) {
        printf("%+f%+fi\n", crealf(output[r]), cimagf(output[r]));
    }
    
    return 0;
}
