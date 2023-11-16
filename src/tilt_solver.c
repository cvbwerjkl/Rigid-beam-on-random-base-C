#include "beam.h"

/*matrices initialization as static for using in all functions without wasting time on reinitialization*/
static gsl_matrix* corMatrix;
static gsl_matrix* springsVector;
static gsl_matrix* springsMatrix;
static gsl_matrix* Vsprings;
static gsl_matrix* nodeCoordinates;
static gsl_matrix* elementNodes;
static gsl_matrix* beamMatrix;
static gsl_matrix* lbeamMatrix; //local  beam element matrix
static gsl_vector* forcesVector;
static gsl_vector* lforcesVector; //local force vector
static gsl_matrix* ranNumVector; //vector of random numbers
static gsl_vector* solution;

static uint32_t elNumGDof = 4; //elements GDof used for local matrixes and vectors creation 
static uint32_t ranNumIndex = 0; //index of random number, used for read data from array of all generated random numbers for all iterations

/*!
 * Формирование корреляционной матрицы через разложение Холецкого и c применениеv корреляционной функции
 * Input:
 * corRadius - радиус корреляцииж.
 * Output:
 * corMatrix - матрица корреляции.
 */
void form_cor_matrix(const uint32_t numberSprings, const double corRadius);


/*!
 * Формирование вектора значений жесткостей пружинок основания с учетом корреляционной матрицы 
 * Input:
 * springStiffnessMean - среднее значение случайной величины жесткости основания;
 * deviation - стандартное отклонение случайной величины жесткости основания;
 * numberSprings - количество пружинок моделирующих основание.
 * Output:
 * springsVector - вектор жесткостей пружинок основания.
 */
void form_springs_vector(const double springStiffnessMean, const double deviation, const uint32_t numberSprings);

/*!
 * Ассемблирования матрицей жесткости стержней с вектором жесткости пружинок и задание граничных условий (закрепляем последний узел по горизонтали) 
 * Input:
 * numberGdof - global degree of freedom;
 * Output:
 * springsMatrix - матрица жесткости пружинок основания.
 */
void form_springs_matrix(const uint32_t numberGdof);

/*!
 * Формирование матрицы жесткости балки 
 * Input:
 * numberElements - количество КЭ балки;
 * ei - изгибная жесткость.
 * Output:
 * beamMatrix - матрица жесткости балки.
 */
void form_beam_el_matrix(const uint32_t numberElements, const double ei);

/*!
 * Формирование вектора сил с заданием граничных условий на последний узел по горизонтали
 * Input:
 * numberElements - количество КЭ балки;
 * pLoad - нагрузка;
 * Output:
 * forcesVector - матрица жесткости балки.
 */
void form_forces_vector(const uint32_t numberElements, const double pLoad);


void tilt_solver(const double deviation, const double springStiffnessMean, const double lBeam, const double corRadius, double* result, uint8_t tests) {
    size_t i;
    uint32_t numbeIter;
    uint32_t numberSprings;
    uint32_t numberElements;
    uint32_t numberNodes;
    uint32_t numberGdof;
    double ei;
    double pLoad;
    double* tilt;
    double* set;
    int s; /*for equtions solving*/
    gsl_permutation* p;

    /* 
     * !!! FOR CRASH TESTS !!! 
     */

    if (tests == 2) { 
        gsl_matrix* matrix;
        gsl_vector* vector;
        gsl_vector* solution;
        matrix = gsl_matrix_alloc(2, 2);
        vector = gsl_vector_alloc(2);
        solution = gsl_vector_alloc(2);
        /*double a = 0;*/
        /*a = gsl_matrix_get(matrix, 2, 0);*/
        /*int c = 10 / 0;*/
        gsl_matrix_set(matrix, 0, 0, 1);
        gsl_matrix_set(matrix, 0, 1, 2);  
        gsl_matrix_set(matrix, 1, 0, 3);
        gsl_matrix_set(matrix, 1, 1, 6);
        gsl_vector_set(vector, 0, 1); 
        gsl_vector_set(vector, 1, 3); 
        
        gsl_permutation* p = gsl_permutation_alloc(2);
        gsl_linalg_LU_decomp (matrix, p, &s);
        gsl_linalg_LU_solve (matrix, p, vector, solution);
        printf ("x = \n");
        gsl_vector_fprintf (stdout, solution, "%g");
        
        gsl_matrix_free(matrix);
        gsl_vector_free(vector);
        gsl_vector_free(solution);
        gsl_permutation_free(p);
    }

    /*
     * !!!! MAIN PART BEGGINING !!!!
     */

    if (tests == 1) { /*parameters for tests*/
        ei = 36 * 100 * 40 * 40 * 40 / 12;
        pLoad = -5;
        numberSprings = 2;
        numbeIter = 1; 
    }
    else { /*parameters for main*/
        ei = 36 * 100 * 40 * 40 * 40 / 12;
        pLoad = -5;
        numberSprings = 100; /*main parameter of beam discretization for FEA*/
        numbeIter = 100;
    }

    numberElements = numberSprings * 2;
    numberNodes = numberElements + 1;
    numberGdof = numberNodes * 2;

    corMatrix = gsl_matrix_alloc(numberSprings, numberSprings);
    springsVector = gsl_matrix_alloc(numberSprings, 1); /*springs stiffness vector*/
    Vsprings = gsl_matrix_alloc(numberSprings, 1); /*temp springs stiffness vector*/
    springsMatrix = gsl_matrix_alloc(numberGdof, numberGdof); /*global stiffness matrix*/
    nodeCoordinates = gsl_matrix_alloc(numberNodes, 1);
    elementNodes = gsl_matrix_alloc(numberElements, 2); /*nodes of elements  !!! переделать в матрицу intов*/
    beamMatrix = gsl_matrix_alloc(numberGdof, numberGdof); 
    lbeamMatrix = gsl_matrix_alloc(4, 4);
    forcesVector = gsl_vector_alloc(numberGdof);  
    lforcesVector = gsl_vector_alloc(4);
    ranNumVector = gsl_matrix_alloc(numberSprings * numbeIter, 1);
    solution = gsl_vector_alloc(numberGdof);
    p = gsl_permutation_alloc(numberGdof);

    tilt = (double*)calloc(numbeIter, sizeof(double));
    set = (double*)calloc(numbeIter, sizeof(double));


    /*nodes coordinates filling in dependance of lBeam and numberNodes*/
    for (i = 0; i < numberNodes; ++i) {
        gsl_matrix_set(nodeCoordinates, i, 0, ((double)lBeam / numberElements) * i);
    }

    /*element nodes filling*/
    for (i = 0; i < numberElements; ++i) {
        gsl_matrix_set(elementNodes, i, 0, i);
        gsl_matrix_set(elementNodes, i, 1, i + 1);
    }

    form_cor_matrix(numberSprings, corRadius);

    if (tests == 1) { /*parameters for tests*/
        gsl_matrix_set(springsVector, 0, 0, 1000);
        gsl_matrix_set(springsVector, 1, 0, 2000);
        //gsl_matrix_set(springsVector, 2, 0, 3);
        //gsl_matrix_set(springsVector, 3, 0, 4);
    }
    else { /*parameters for main*/
        /*Random numbers vector creation start*/
        const gsl_rng_type* T;
        gsl_rng* r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        
        gsl_rng_set(r, time(NULL)); /*for different random numbers in dependence of the current time*/
        
        /*random numers generation for all iterations*/
        for (i = 0; i < numberSprings * numbeIter; ++i) {
            gsl_matrix_set(ranNumVector, i, 0, gsl_ran_gaussian(r, deviation));
        }
        gsl_rng_free(r);
        
        FILE * dd = fopen("testranNum.txt", "wb"); /*saving random numbers for checking*/
        gsl_matrix_fprintf(dd, ranNumVector, "%f");
        fclose (dd);
        /*Random numbers vector creation end*/
    }

    /*matrix generation for all iteractions*/
    form_beam_el_matrix(numberElements, ei);
    form_forces_vector(numberElements, pLoad);

    for (i = 0; i < numbeIter; ++i) {
        /*matrix generation with random*/
        if (tests != 1) {
            form_springs_vector(springStiffnessMean, deviation, numberSprings);
        }
        form_springs_matrix(numberGdof);

        /*solving*/
        gsl_linalg_LU_decomp (springsMatrix, p, &s);
        gsl_linalg_LU_solve (springsMatrix, p, forcesVector, solution);

        tilt[i] = fabs(gsl_vector_get(solution, 0) - gsl_vector_get(solution,numberGdof - 2)) / lBeam;   
        set[i] = fabs(gsl_vector_get(solution, 0) + gsl_vector_get(solution,numberGdof - 2)) / 2;   
    }

    /*
     * Statistical processing of results
     */
    result[0] = gsl_stats_mean(tilt, 1, numbeIter);
    result[1] = gsl_stats_max(tilt, 1, numbeIter);
    result[2] = gsl_stats_min(tilt, 1, numbeIter);
    if (tests != 1) {
        result[3] = gsl_stats_sd_m(tilt, 1, numbeIter, result[0]);
        result[4] = result[0] * gsl_cdf_tdist_Pinv(1 - (1 - 0.95) / 2, numbeIter - 1); /*tilt with 0.95 confidence*/
        result[5] = result[0] * gsl_cdf_tdist_Pinv(1 - (1 - 0.99) / 2, numbeIter - 1); /*tilt with 0.99 confidence*/
        result[6] = result[4] * (result[3] / (sqrt((double)numbeIter))); /*error for 0.95 confidence*/
        result[7] = result[5] * (result[3] / (sqrt((double)numbeIter))); /*error for 0.95 confidence*/
    }
    result[8] = gsl_stats_mean(set, 1, numbeIter);


    if (tests == 1) {
        FILE * cM = fopen("testCorMatrix.txt", "wb");
        gsl_matrix_fprintf(cM, corMatrix, "%f");
        fclose (cM);
        FILE * sV = fopen("testSpringsVector.txt", "wb");
        gsl_matrix_fprintf(sV, springsVector, "%f");
        fclose (sV);
        FILE * bM = fopen("testBeamMatrix.txt", "wb");
        gsl_matrix_fprintf(bM, beamMatrix, "%f");
        fclose (bM);
        FILE * fV = fopen("testforcesVector.txt", "wb");
        gsl_vector_fprintf(fV, forcesVector, "%f");
        fclose (fV);
        FILE * sM = fopen("testSpringsMatrix.txt", "wb");
        gsl_matrix_fprintf(sM, springsMatrix, "%f");
        fclose (sM);
        FILE * sO = fopen("testsolutionVector.txt", "wb");
        gsl_vector_fprintf(sO, solution, "%f");
        fclose (sO);
    }
    
    gsl_matrix_free(corMatrix);
    gsl_matrix_free(springsVector);
    gsl_matrix_free(Vsprings);
    gsl_matrix_free(springsMatrix);
    gsl_matrix_free(nodeCoordinates);
    gsl_matrix_free(elementNodes);
    gsl_matrix_free(beamMatrix);
    gsl_matrix_free(lbeamMatrix);
    gsl_vector_free(forcesVector);
    gsl_vector_free(lforcesVector);
    gsl_matrix_free(ranNumVector);
    gsl_vector_free(solution);
    gsl_permutation_free(p);
    free(tilt);
    free(set);

    return;
}

void form_cor_matrix(const uint32_t numberSprings, const double corRadius) {
    size_t i;
    size_t j;

    for (i = 0; i < numberSprings; ++i) {
        for (j = 0; j < numberSprings; ++j) {
          gsl_matrix_set(corMatrix, i, j, exp(-2 * fabs(gsl_matrix_get(nodeCoordinates, j * 2, 0) - gsl_matrix_get(nodeCoordinates, i * 2, 0)) / corRadius));
        }
    }

    /*Cholesky decomposition A = L * Ltran*/
    gsl_linalg_cholesky_decomp1(corMatrix);

    /*clening unwanted nonzero values from L matrix*/
    for (i = 0; i < numberSprings - 1; i++) {
      for (j = i + 1; j < numberSprings; j++) {
       gsl_matrix_set(corMatrix, i, j, 0);
      }
    }

    return;
}

void form_springs_vector(const double springStiffnessMean, const double deviation, const uint32_t numberSprings) {
    size_t i;
    double a; /*normal distribution parameter mean with elements amount consideration*/
    double b; /*normal distribution parameter standard deviation with elements amount consideration*/

    a = springStiffnessMean / (double)numberSprings;
    b = (springStiffnessMean * deviation) / (double)numberSprings;

    
    for (i = 0; i < numberSprings; ++i) {
        gsl_matrix_set(Vsprings, i, 0, gsl_matrix_get(ranNumVector, ranNumIndex, 0)); 
        ++ranNumIndex;
    }
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, corMatrix, Vsprings, 0.0, springsVector); /*matrices multiplication */
    gsl_matrix_scale(springsVector, b);
    gsl_matrix_add_constant(springsVector, a);

    return;
}

void form_springs_matrix(const uint32_t numberGdof) {
    size_t i;
    size_t j;
    int k = 0;

    for (i = 0; i < numberGdof; i++) {
      for (j = 0; j < numberGdof; j++) {
       gsl_matrix_set(springsMatrix, i, j, gsl_matrix_get(beamMatrix, i, j));
      }
    }
    for (i = 2; i < numberGdof; i += 4) {
        gsl_matrix_set(springsMatrix, i, i, gsl_matrix_get(springsMatrix, i, i) + fabs(gsl_matrix_get(springsVector, k, 0)));
        ++k;
    }

    return;
}

void form_beam_el_matrix(const uint32_t numberElements, const double ei) {
    size_t i;
    size_t j;
    size_t k;
    size_t node0;
    size_t node1;
    uint32_t elementDof[elNumGDof];
    double elemLength;

    for (i = 0; i < numberElements; ++i) {
        node0 = gsl_matrix_get(elementNodes, i, 0);
        node1 = gsl_matrix_get(elementNodes, i, 1);
        elementDof[0] = 2 * node0;
        elementDof[1] = 2 * node1 - 1;
        elementDof[2] = 2 * node1;
        elementDof[3] = 2 * node1 + 1;
        elemLength = gsl_matrix_get(nodeCoordinates, node1, 0) - gsl_matrix_get(nodeCoordinates, node0, 0);

        /*local beam elements matrix filling*/
        gsl_matrix_set(lbeamMatrix, 0, 0, 12);
        gsl_matrix_set(lbeamMatrix, 0, 1, 6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 0, 2, -12);
        gsl_matrix_set(lbeamMatrix, 0, 3, 6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 1, 0, 6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 1, 1, 4 * elemLength * elemLength);
        gsl_matrix_set(lbeamMatrix, 1, 2, -6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 1, 3, 2 * elemLength * elemLength);
        gsl_matrix_set(lbeamMatrix, 2, 0, -12);
        gsl_matrix_set(lbeamMatrix, 2, 1, -6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 2, 2, 12);
        gsl_matrix_set(lbeamMatrix, 2, 3, -6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 3, 0, 6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 3, 1, 2 * elemLength * elemLength);
        gsl_matrix_set(lbeamMatrix, 3, 2, -6 * elemLength);
        gsl_matrix_set(lbeamMatrix, 3, 3, 4 * elemLength * elemLength);

        gsl_matrix_scale(lbeamMatrix, (double)ei / (elemLength * elemLength * elemLength));

        /*assembling global beam elements matrix*/
        for (j = 0; j < elNumGDof; ++j) { 
            for (k = 0; k < elNumGDof; ++k) {
                gsl_matrix_set(beamMatrix, elementDof[j], elementDof[k], gsl_matrix_get(beamMatrix, elementDof[j], elementDof[k]) + gsl_matrix_get(lbeamMatrix, j, k));
            }
        }
    }

    return; 
}

void form_forces_vector(const uint32_t numberElements, const double pLoad) {
    size_t i;
    size_t j;
    size_t node0;
    size_t node1;
    uint32_t elementDof[4];
    double elemLength;

    for (i = 0; i < numberElements; ++i) {
        node0 = gsl_matrix_get(elementNodes, i, 0);
        node1 = gsl_matrix_get(elementNodes, i, 1);
        elementDof[0] = 2 * node0;
        elementDof[1] = 2 * node1- 1;
        elementDof[2] = 2 * node1;
        elementDof[3] = 2 * node1 + 1;
        elemLength = gsl_matrix_get(nodeCoordinates, node1, 0) - gsl_matrix_get(nodeCoordinates, node0, 0);
    
        /*local force vector filling*/
        gsl_vector_set(lforcesVector, 0, (pLoad * elemLength) / 2);
        gsl_vector_set(lforcesVector, 1, (pLoad * elemLength * elemLength) / 12);
        gsl_vector_set(lforcesVector, 2, (pLoad * elemLength) / 2);
        gsl_vector_set(lforcesVector, 3, (-pLoad * elemLength * elemLength) / 12);

        /*assembling global forces vector*/
        for (j = 0; j < elNumGDof; ++j) {
            gsl_vector_set(forcesVector, elementDof[j], gsl_vector_get(forcesVector, elementDof[j]) + gsl_vector_get(lforcesVector, j));
        }
        
    }

    return;
}
