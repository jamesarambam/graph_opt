#include <stdio.h>
#include <stdlib.h>
#include "rcdcop.h"
#include <time.h>

// -------------------------------------------- //

double gateWay(int nT, int dT, double *t1, double (*t2a)[nT], double (*t2b)[nT], double (*pot)[nT][dT][dT],
               double (*theta)[nT][dT][dT], double (*thetaH)[nT][dT][dT], double (*Pij)[nT][dT][dT], double (*Pstar_ij)[nT][dT][dT], int totIter, int mR) {

    printf("# ============= C ================ #\n");

    clock_t begin = clock();
    // ------------------ Variables ------------------ //
    int N = nT;
    int D[N];
    int nbr_l[N][N];
    int nbr_q[N][N];
    struct node Nodes[N];

    mDomain = dT;
    tNodes = nT;
    maxRuntime = mR;

    // ------------------------- Type Casting ----------------- //

    for (int i = 0; i < N; ++i) {
        D[i] = (int) t1[i];
    }
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            nbr_l[j][i] = (int) t2a[j][i];
        }
    }
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            nbr_q[j][i] = (int) t2b[j][i];
        }
    }

    // ----------------- Nodes -------------- //
    for (int i = 0; i < N; ++i) {
        Nodes[i].domain = D[i];
    }

    // ----------------- Edges -------------- //
    int totE, qpLen, lpLen;
    totE = N * (N - 1) / 2;
    struct edge QP[totE];
    for (int i = 0; i < totE; ++i) {
        QP[i].i = -1;
        QP[i].j = -1;
    }
    qpLen = getE(N, nbr_q, QP);

    struct edge LP[totE];
    for (int i = 0; i < totE; ++i) {
        LP[i].i = -1;
        LP[i].j = -1;
    }
    lpLen = getE(N, nbr_l, LP);

    // --------------------------------------- //

    double Pstar[nT][dT];
    double PiFi[nT][dT];
    double Li[nT];
    double Lijxi[nT][nT][dT];

    for (int i = 0; i < nT; ++i) {
        for (int j = 0; j < nT; ++j) {
            for (int xi = 0; xi < Nodes[i].domain; ++xi) {
                Lijxi[i][j][xi] = 0;
            }
        }
    }


    // --------------------------------------- //

    FILE *f = fopen("file.log", "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    double intObj, total_time = 0;
    initialize(nT, dT, Pstar, Pstar_ij, Nodes, lpLen, LP);

    iter_count = 1;

    while (iter_count <= totIter) {
//    while (true) {

        clock_t t1 = clock();
        printf("# ------------------------- #\n");
        printf("RCDCOP Iteration : %d\n", iter_count);
        fprintf(f,"# ---------------------- #\n");
        fprintf(f,"Iteration : %d\n", iter_count);

        updatePij(Pij, Pstar_ij);
        checkThreshold_P(Nodes, Pstar);
        checkThreshold_Pij(Pij);
        updatePiFi(nT, dT, Nodes, Pstar, PiFi, thetaH, nbr_q);
        solveBCD(PiFi, Li, Lijxi, Nodes, nbr_l, Pij, lpLen, LP, thetaH, f);
        computePstar(Nodes, PiFi, Li, Lijxi, nbr_l, Pstar);
        computePstar_ij(Pij, Pstar_ij, thetaH, Lijxi, lpLen, LP, Nodes);
        printf(" Total BCD Iteration : %d\n", inner_iteration);
        printf(" BCD Runtime : %f\n", bcdRuntime);
        fprintf(f," Total BCD Iteration : %d\n", inner_iteration);
        fprintf(f," BCD Runtime : %f\n", bcdRuntime);
        double cObj = computeObjective(pot, Pstar, Pstar_ij, qpLen, QP, lpLen, LP, Nodes);
        int sol[tNodes];
        primalExtract(Pstar, sol, nbr_l, nbr_q, Nodes, theta);
        intObj = computeIntObjective(pot, qpLen, QP, lpLen, LP, sol);

        clock_t t2 = clock();
        double itr_runtime = (double)(t2 - t1) / CLOCKS_PER_SEC;
        total_time += itr_runtime;

        printf(" Continuous Objective : %f\n", cObj);
        printf(" Integer Objective : %f\n", intObj);
        printf(" Iteration Runtime : %f\n", itr_runtime);
        printf(" Total Runtime : %f\n", total_time);

        fprintf(f, " Continuous Objective : %f\n", cObj);
        fprintf(f, " Integer Objective : %f\n", intObj);
        fprintf(f, " Iteration Runtime : %f\n", itr_runtime);
        fprintf(f, " Total Runtime : %f\n", total_time);

        iter_count += 1;

        if (total_time > maxRuntime){
            break;
        }

    }

    // -------------------------------------- //

    double Obj = computeObjective(pot, Pstar, Pstar_ij, qpLen, QP, lpLen, LP, Nodes);
    printf("# ================== Final Result ===================== #\n");
    printf("Total Outer Iteration : %d\n", iter_count-1);
    printf("Continuous Objective : %f\n", Obj);
    printf("Integer Objective : %f\n", intObj);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total Execution Time in C : %f Seconds\n", time_spent);

    fprintf(f,"# ================== Final Result ===================== #\n");
    fprintf(f, "Total Outer Iteration : %d\n", iter_count-1);
    fprintf(f, "Continuous Objective : %f\n", Obj);
    fprintf(f, "Integer Objective : %f\n", intObj);
    fprintf(f,"Total Execution Time in C : %f Seconds\n", time_spent);

    fclose(f);


    return 0.0;
}


double test(int *t1) {


    printf("============= C ================\n");

    printf("%d\n", t1[0]);
    printf("%d\n", t1[1]);
    printf("%d\n", t1[2]);

    return 0.0;
}
