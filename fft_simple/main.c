#include <stdio.h>
#include <math.h>
#include <complex.h>


void generate_twiddle_factor_lookup_table(int n,
                                          complex double twiddle_lut[]) {
    int i;
    for (i = 0; i < n; i++) {
        twiddle_lut[i] = cexp(((-I * 2 * M_PI) / n) * i);
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
    printf("bit-reversing...\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < (int)log2(n); j++) {
            bits[j] = i >> j & 0x1;
        }
        reversed = 0;
        for (j = 0; j < (int)log2(n); j++) {
            reversed = reversed | (bits[(int)log2(n) - 1 - j] << j);
        }
        printf("[%2d] -> [%2d]\n", i, reversed);
        y[reversed] = x[i];
    }
}

void butterfly_block(complex double x0,
                     complex double x1,
                     complex double wk,
                     complex double *y0,
                     complex double *y1) {
    *y0 = x0 + wk * x1;
    *y1 = x0 - wk * x1;
}


int main(int argc, const char * argv[]) {
    
    // define transform size (n needs to be a power of 2)
    int n = 16;
    
    // set input
    complex double input[n];
    input[ 0] = +0.0 + 0.0*I;
    input[ 1] = +0.0 + 0.0*I;
    input[ 2] = +1.0 + 0.0*I;
    input[ 3] = +1.0 + 0.0*I;
    input[ 4] = +0.0 + 0.0*I;
    input[ 5] = +0.0 + 0.0*I;
    input[ 6] = +1.0 + 0.0*I;
    input[ 7] = +1.0 + 0.0*I;
    input[ 8] = +0.0 + 0.0*I;
    input[ 9] = +0.0 + 0.0*I;
    input[10] = +1.0 + 0.0*I;
    input[11] = +1.0 + 0.0*I;
    input[12] = +0.0 + 0.0*I;
    input[13] = +0.0 + 0.0*I;
    input[14] = +1.0 + 0.0*I;
    input[15] = +1.0 + 0.0*I;
    
    // generate twiddle factor look up table
    complex double twiddle_lut[n];
    printf("generating twiddle factors...\n");
    generate_twiddle_factor_lookup_table(n, twiddle_lut);
    
    // generate butterfly blocks and all the interconnections
    complex double stage_wires[(int)log2(n)][n];
    bit_reverse_block(n, input, stage_wires[0]);
    int stage, group_offset, butterfly_offset;
    printf("generating buttefly structure...\n");
    
    for (stage = 1;
         stage <= (int)log2(n);
         stage++) {
        
        for (group_offset = 0;
             group_offset <= n-1;
             group_offset = group_offset + (int)pow(2.0, stage)) {
            
            for (butterfly_offset = 0;
                 butterfly_offset <= (int)pow(2.0, stage)/2-1;
                 butterfly_offset++) {
                
                butterfly_block(stage_wires[stage-1][group_offset+butterfly_offset],
                                stage_wires[stage-1][group_offset+butterfly_offset+(int)pow(2.0, stage)/2],
                                twiddle_lut[(int)pow(2.0, (int)log2(n)-stage)*butterfly_offset],
                                &stage_wires[stage][group_offset+butterfly_offset],
                                &stage_wires[stage][group_offset+butterfly_offset+(int)pow(2.0, stage)/2]);
                
                printf("stage: %d,\tgroup_offset: %d,\tbutterfly_offset: %d,\ttwiddle_factor: %d\n",
                       stage,
                       group_offset,
                       butterfly_offset,
                       (int)pow(2.0, (int)log2(n)-stage)*butterfly_offset);
            }
        }
    }
    
    // connect output to the last stage
    complex double output[n];
    int i;
    for (i = 0; i < n; i++) {
        output[i] = stage_wires[(int)log2(n)][i];
    }
    
    
    printf("input vector:\n");
    for (i = 0; i < n; i++) {
        printf("[%2d] %+f%+fi\n",
               i,
               crealf(input[i]),
               cimagf(input[i]));
    }
    
    printf("bit-reversed input vector:\n");
    for (i = 0; i < n; i++) {
        printf("[%2d] %+f%+fi\n",
               i,
               crealf(stage_wires[0][i]),
               cimagf(stage_wires[0][i]));
    }
    
    printf("output vector:\n");
    for (i = 0; i < n; i++) {
        printf("[%2d] %+f%+fi\n",
               i, crealf(output[i]),
               cimagf(output[i]));
    }
    
    return 0;
}
