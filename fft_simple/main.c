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
#include <string.h>


void generate_twiddle_factor_lookup_table(int n,
                                          complex double twiddle_lut[]) {
    int i;
    for (i = 0; i < n; i++) {
        twiddle_lut[i] = cexp(((-I * 2 * M_PI) / n) * i);
        // TODO: only using first half of the table but when shortened, we get SEGFAULTS...
        if (i > n/2-1) {
            twiddle_lut[i] = 0.0+0.0*I;
        }
        printf("%+f%+fi\n",
               crealf(twiddle_lut[i]),
               cimagf(twiddle_lut[i]));
    }
}

void bit_reverse_block(int n,
                       complex double x[],
                       complex double y[]) {
    int i, j;
    int reversed;
    int bits[(int)log2(n)];
    //printf("bit-reversing...\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < (int)log2(n); j++) {
            bits[j] = i >> j & 0x1;
        }
        reversed = 0;
        for (j = 0; j < (int)log2(n); j++) {
            reversed = reversed | (bits[(int)log2(n) - 1 - j] << j);
        }
        //printf("[%2d] -> [%2d]\n", i, reversed);
        y[reversed] = x[i];
    }
}

void butterfly_block(complex double x0,
                     complex double x1,
                     complex double wk,
                     complex double *y0,
                     complex double *y1) {
    complex double wkx1;
    wkx1 = wk * x1;
    *y0 = x0 + wkx1;
    *y1 = x0 - wkx1;
}


void dit_radix_2_fft(int n, complex double input[], complex double output[]) {
    
    // generate twiddle factor look up table
    complex double twiddle_lut[n];
    //printf("generating twiddle factors...\n");
    generate_twiddle_factor_lookup_table(n, twiddle_lut);
    
    // generate butterfly blocks and all the interconnections
    complex double stage_wires[(int)log2(n)][n];
    bit_reverse_block(n, input, stage_wires[0]);
    int stage, group_offset, butterfly_offset;
    //printf("generating buttefly structure...\n");
    
    for (stage = 1;
         stage <= (int)log2(n);
         stage++) {
        
        for (group_offset = 0;
             group_offset <= n-1;
             group_offset = group_offset + (int)pow(2.0, stage)) {
            
            for (butterfly_offset = 0;
                 butterfly_offset <= (int)pow(2.0, stage)/2-1;
                 butterfly_offset++) {
                
                int x0_i, x1_i, wk_i, y0_i, y1_i;
                
                x0_i = group_offset+butterfly_offset;
                x1_i = group_offset+butterfly_offset+(int)pow(2.0, stage)/2;
                wk_i = (int)pow(2.0, (int)log2(n)-stage)*butterfly_offset;
                y0_i = group_offset+butterfly_offset;
                y1_i = group_offset+butterfly_offset+(int)pow(2.0, stage)/2;
                
                printf("stage: %d,\tgroup_offset: %d,\tbutterfly_offset: %d,\ttwiddle_factor: %d\n",
                       stage,
                       group_offset,
                       butterfly_offset,
                       wk_i);
                
                butterfly_block(stage_wires[stage-1][x0_i],
                                stage_wires[stage-1][x1_i],
                                twiddle_lut[wk_i],
                                &stage_wires[stage][y0_i],
                                &stage_wires[stage][y1_i]);
                

            }
        }
    }
    
    // connect output to the last stage
    memcpy(output, stage_wires[(int)log2(n)], sizeof(complex double)*n);
}

void check_result(int n, double eps, double complex expected[], double complex actual[]) {
    int i;
    for (i = 0; i < n; i++) {
        if (fabs(cabs(expected[i])-cabs(actual[i])) > eps) {
            printf("failed\n");
            return;
        }
    }
    printf("passed\n");
}



int main(int argc, const char * argv[]) {
    
    // define transform size (n needs to be a power of 2)
    int n = 16;
    
    // set input
    complex double input[n];
    input[ 0] =+0.000000e+00;
    input[ 1] =+0.95106e+00;
    input[ 2] =+0.58779e+00;
    input[ 3] =-0.58779e+00;
    input[ 4] =-0.95106e+00;
    input[ 5] =-1.1331e-15;
    input[ 6] =+0.95106e+00;
    input[ 7] =+0.58779e+00;
    input[ 8] =-0.58779e+00;
    input[ 9] =-0.95106e+00;
    input[10] =-2.26620e-15;
    input[11] =+0.95106e+00;
    input[12] =+0.58779e+00;
    input[13] =-0.58779e+00;
    input[14] =-0.95106e+00;
    input[15] =-4.28750e-15;
    
    // prepare output array
    complex double output[n];
    
    // build fft
    dit_radix_2_fft(n, input, output);
    
    printf("input vector:\n");
    int i;
    for (i = 0; i < n; i++) {
        printf("%+f%+fi\n",
               crealf(input[i]),
               cimagf(input[i]));
    }
    
    printf("output vector:\n");
    for (i = 0; i < n; i++) {
        printf("%+f%+fi\n",
               crealf(output[i]),
               cimagf(output[i]));
    }
    
    complex double expected[16];
    expected[ 0] =-2.664535e-15+0.000000e+00*I;
    expected[ 1] =+5.887077e-02-2.959633e-01*I;
    expected[ 2] =+3.498683e-01-8.446568e-01*I;
    expected[ 3] =+3.984873e+00-5.963785e+00*I;
    expected[ 4] =-1.538841e+00+1.538841e+00*I;
    expected[ 5] =-9.505633e-01+6.351460e-01*I;
    expected[ 6] =-7.988962e-01+3.309136e-01*I;
    expected[ 7] =-7.420404e-01+1.476010e-01*I;
    expected[ 8] =-7.265425e-01+0.000000e+00*I;
    expected[ 9] =-7.420404e-01-1.476010e-01*I;
    expected[10] =-7.988962e-01-3.309136e-01*I;
    expected[11] =-9.505633e-01-6.351460e-01*I;
    expected[12] =-1.538841e+00-1.538841e+00*I;
    expected[13] =+3.984873e+00+5.963785e+00*I;
    expected[14] =+3.498683e-01+8.446568e-01*I;
    expected[15] =+5.887077e-02+2.959633e-01*I;
    
    printf("result verification:\n");
    check_result(n, 1.0e-4, expected, output);
    
    return 0;
}
