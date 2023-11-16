#include "beam.h"

int main(int argc, char **argv) {
    double deviation;
    double springStiffnessMean;
    double lBeam;
    double corRadius;
    double result[9] = {0};

    uint8_t tests;
    char param;

    param = getopt(argc, argv, "tc");
    if (param == 255) {
        tests = 0;
    }
    else if (param == 't') {
        tests = 1;
    }
    else if (param == 'c') {
        tests = 2;
    }
    else {
        fprintf(stderr, "Unknown option `-%c'.\n", param);
        return 1;
    }

    printf("Realisation type %d\n", tests);

    if (tests == 1) { /*paramteters for tests*/
        deviation = 0.15;
        springStiffnessMean = 1000.0;
        lBeam = 25.0;
        corRadius = 1.0;
    }
    else { /*paramteters for main*/
        deviation = 0.15;
        springStiffnessMean = 1000.0;
        lBeam = 25.0;
        corRadius = 1.0;
    }

    tilt_solver(deviation, springStiffnessMean, lBeam, corRadius, result, tests);
        
    printf("Avg tilt %f\n", result[0]);
    printf("Max tilt %f\n", result[1]);
    printf("Min tilt %f\n", result[2]);
    printf("Tilt standart deviation %f\n", result[3]);
    printf("Tilt 0.95 confidence %f\n", result[4]);
    printf("Tilt 0.99 confidence %f\n", result[5]);
    printf("Error 0.95 confidence %.10f\n", result[6]);
    printf("Error 0.99 confidence %.10f\n", result[7]);
    printf("Average settlement %f\n", result[8]);
    
    return 0;
}
