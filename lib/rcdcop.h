//
// Created by James Arambam on 22/8/16.
//

#ifndef CLIONPROJECTS_RCDCOP_H
#define CLIONPROJECTS_RCDCOP_H

#include <stdlib.h>
#include "gslBrent.h"
#include "gslNewton.h"
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>

// ----------------------------- //
double zeroThreshold = 1e-10;
double P_Threshold = 1e-5;
double conVioThres = 0.001;
int mDomain, maxRuntime;
int tNodes;
int inner_iteration, totRuntime, iter_count;
double bcdRuntime;
int bcdConv = 20, bcdConflag1 = 0, bcdConflag2 = 0;
double oldEst1 = -1, oldEst2 = -1;
// ---------------------------- //

//void write(FILE *f, )

void halt() {
    exit(0);
}

double roundDec(double x, unsigned int digits) {
    double fac = pow(10, digits);
    double t1 = x * fac;
    int t2 = (int) t1;
    return (double)t2 / fac;
}

double max1D(int size, double *arr) {

    double max = -100.0;

    for (int i = 0; i < size; ++i) {
        if (max <= arr[i]) {
            max = arr[i];
        }
    }
    return max;
}

struct edge {
    int i;
    int j;
};

struct node {
    int domain;
};

int getE(int n, int (*nq)[n], struct edge *E) {

    int eC = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (nq[i][j] != -1 && nq[i][j] > i) {
                E[eC].i = i;
                E[eC].j = nq[i][j];
                eC = eC + 1;
            }
        }
    }
    return eC;
}

void initialize(int n, int d, double (*P)[d], double (*Pij)[n][d][d], struct node *Nodes, int lenLP, struct edge *LP) {


    // --------------- Pi ----------------- //
    for (int i = 0; i < n; ++i) {
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            P[i][xi] = 1 / (double) Nodes[i].domain;
        }
    }

    // --------------- Pij ---------------- //
    for (int e = 0; e < lenLP; ++e) {
        int i = LP[e].i;
        int j = LP[e].j;
        if (i != -1 && j != -1) {
            for (int xi = 0; xi < Nodes[i].domain; ++xi) {
                for (int xj = 0; xj < Nodes[j].domain; ++xj) {
                    double pij = P[i][xi] * P[j][xj];
                    Pij[i][j][xi][xj] = pij;
                    Pij[j][i][xj][xi] = pij;
                }
            }
        }
    }
}

double f(int i, int xi, int nT, int d, double (*P)[d], double (*thetaH)[nT][d][d], int *nbr_q, struct node *Nodes) {

    double val = 0;
    int j;
    for (int n = 0; n < nT; ++n) {
        j = nbr_q[n];
        if (j != -1) {
            for (int xj = 0; xj < Nodes[j].domain; ++xj) {
                val = val + (thetaH[i][j][xi][xj] * P[j][xj]);
            }
        }
    }
    return val;
}

void updatePiFi(int n, int d, struct node *Nodes, double (*P)[d], double (*PiFi)[d], double (*thetaH)[n][d][d],
                int (*nbr_q)[n]) {

    for (int i = 0; i < n; ++i) {
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            PiFi[i][xi] = P[i][xi] * f(i, xi, n, d, P, thetaH, nbr_q[i], Nodes);
        }
    }
}

void init_Li_Lijxi(int n, int d, int (*nbr_l)[n], struct node *Nodes, double (*Lijxi)[n][d], double *Li) {

    int nb_j;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            nb_j = nbr_l[i][j];
            if (nb_j != -1) {
                for (int xi = 0; xi < Nodes[i].domain; ++xi) {
                    Lijxi[i][nb_j][xi] = 1;
                }
            }

        }

    }

    double m;
    for (int i = 0; i < n; ++i) {
        double tmp[Nodes[i].domain];
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            int tot = 0;
            for (int j = 0; j < n; ++j) {
                nb_j = nbr_l[i][j];
                tot = tot + Lijxi[i][nb_j][xi];
            }
            tmp[xi] = tot;
        }
        m = max1D(Nodes[i].domain, tmp);
        Li[i] = m + 2;
    }

}

double computeSum_lijxi(int i, int xi, double (*Lijxi)[tNodes][mDomain], int (*nbr_l)[tNodes]) {

    double tot = 0;
    int nb_j;
    for (int j = 0; j < tNodes; ++j) {
        nb_j = nbr_l[i][j];
        if (nb_j != -1) {
            tot = tot + Lijxi[i][nb_j][xi];
            // printf("%d %d %d %f\n", i, nb_j, xi, Lijxi[i][nb_j][xi]);

        }
    }
    return tot;
}

void getLiParam(int n, int d, int i, double *fP, double (*PiFi)[d], double (*Lijxi)[n][d], struct node *Nodes,
                int (*nbr_l)[n]) {

    int t_ind = 0;
    for (int xi = 0; xi < Nodes[i].domain; ++xi) {
        fP[t_ind] = PiFi[i][xi];
        fP[t_ind + 1] = -1 * computeSum_lijxi(i, xi, Lijxi, nbr_l);
//        printf("%d, %f\n", t_ind, PiFi[i][xi]);
//        printf("%d, %f\n", t_ind+1, -1 * computeSum_lijxi(i, xi, Lijxi, nbr_l));

        t_ind = t_ind + 2;
    }
    fP[2 * Nodes[i].domain] = -1;

}

int getBounds_Li(int n, int d, int i, int pLen, double *fP, double (*Lijxi)[n][d], double *bounds, struct node *Nodes,
                 int (*nbr_l)[n]) {

    struct func_Params p;
    double lb_tmp, eps = 0.15, x, lb, ub, tmp[Nodes[i].domain];
    int count = 0, it = 0, mIter = 500;
    bool flag;

    flag = true;
    p.d = fP;
    param_len = pLen;
    // --------------------------------- //
    for (int xi = 0; xi < Nodes[i].domain; ++xi) {
        tmp[xi] = computeSum_lijxi(i, xi, Lijxi, nbr_l);
    }
    lb_tmp = max1D(Nodes[i].domain, tmp);
    x = lb_tmp + eps;
    ub = x;
    // -------------------------------- //
    while (true) {
        if (fabs(ratFun(x, (void *) &p)) <= zeroThreshold) {
            bounds[0] = x;
            bounds[1] = x;
            return 1;
        }
        if (ratFun(x, (void *) &p) > 0) {
            if (flag) {
                lb = x;
                x = lb + eps;
                ub = x;
            } else {
                if (ratFun(x, (void *) &p) > 0 && ratFun(ub, (void *) &p) < 0) {
                    bounds[0] = x;
                    bounds[1] = ub;
                    return 1;
                } else {
                    printf("Bounds Not Found, Exiting ! for Li %d\n", i);
                    halt();
                }
            }
        } else {
            flag = false;
            count = count + 1;
            x = x - eps / pow(2, count);
        }
        it = it + 1;
        if (it > mIter) {
            printf("Unable to find bounds for Li %d\n", i);
            return 0;
        }
    }
}

void updateLambda_i(int n, int d, double (*PiFi)[d], double *Li, double (*Lijxi)[n][d], struct node *Nodes,
                    int (*nbr_l)[n]) {

    double root;
    double bounds[2], lb, ub;

    for (int i = 0; i < n; ++i) {
        double fP[2 * Nodes[i].domain + 1];
        bounds[0] = 0;
        bounds[1] = 0;
        getLiParam(n, d, i, fP, PiFi, Lijxi, Nodes, nbr_l);
        int pLen = 2 * Nodes[i].domain + 1;
        getBounds_Li(n, d, i, pLen, fP, Lijxi, bounds, Nodes, nbr_l);
        lb = bounds[0];
        ub = bounds[1];
        if (lb == ub) {
            Li[i] = lb;
        } else {
            root = getRoot(fP, pLen, lb, ub);
            Li[i] = root;
            // printf("New Root 1 | ");
        }
    }
}

void computePstar(struct node *Nodes, double (*PiFi)[mDomain], double *Li, double (*Lijxi)[tNodes][mDomain],
                  int (*nbr_l)[tNodes], double (*Pstar)[mDomain]) {

    double Den;
    for (int i = 0; i < tNodes; ++i) {
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            Den = Li[i] - computeSum_lijxi(i, xi, Lijxi, nbr_l);
            Pstar[i][xi] = PiFi[i][xi] / Den;
        }
    }
}

double computeViolation_Pi(struct node *Nodes, double (*PiFi)[mDomain], double *Li, double (*Lijxi)[tNodes][mDomain],
                           int (*nbr_l)[tNodes]) {

    // ------- Compute Pstar ------- /
    double Pstar[tNodes][mDomain];

    computePstar(Nodes, PiFi, Li, Lijxi, nbr_l, Pstar);


//     printf("-------------\n");
//     for (int i = 0; i < tNodes; ++i) {
//         for (int xi = 0; xi < Nodes[i].domain; ++xi) {
//             printf("%d %d %.15f\n", i, xi, Pstar[i][xi]);
//         }
//     }
//     printf("-------------\n");

    // ---- Compute Violation ----- //
    double temp[tNodes];
    for (int i = 0; i < tNodes; ++i) {
        double sum = 0;
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            sum += Pstar[i][xi];
        }
        sum -= 1;
        // temp[i] = fabs(sum);
        temp[i] = sum;
    }


//    printf("-------------\n");
//    for (int j = 0; j < tNodes; ++j) {
//        printf("%d : %.16f\n", j, temp[j]);
//    }
//    printf("-------------\n");

    return max1D(tNodes, temp);
}

double computeSum_lijxi_k(int i, int j, int xi, double (*Lijxi)[tNodes][mDomain], int (*nbr_l)[tNodes]) {

    double tot = 0;
    int nb_k;
    for (int k = 0; k < tNodes; ++k) {
        nb_k = nbr_l[i][k];
        if (nb_k != -1 && nb_k != j) {
            tot = tot + Lijxi[i][nb_k][xi];
        }
    }
    return tot;
}

void getLijxiParam(int i, int j, int xi, double *fP, double (*PiFi)[mDomain], double (*Pstar_ij)[tNodes][mDomain][mDomain],
              double *Li, double (*Lijxi_old)[tNodes][mDomain], struct node *Nodes, int (*nbr_l)[tNodes],
              double (*thetaH)[tNodes][mDomain][mDomain]) {

    int t_ind = 0;
    int fp_size = 2 * Nodes[j].domain + 2;
    for (int xj = 0; xj < Nodes[j].domain; ++xj) {

//        printf("%d %d %d %d %f\n", i, j, xi, xj, thetaH[i][j][xi][xj]);
//        printf("%d %d %d %d %f\n", i, j, xi, xj, Pstar_ij[i][j][xi][xj]);
        fP[t_ind] = thetaH[i][j][xi][xj] * Pstar_ij[i][j][xi][xj];
        fP[t_ind + 1] = Lijxi_old[j][i][xj];
        t_ind += 2;
    }
    fP[fp_size - 2] = PiFi[i][xi];
    fP[fp_size - 1] = -1 * Li[i] + computeSum_lijxi_k(i, j, xi, Lijxi_old, nbr_l);
}

int computeBounds_lijxi(int i, int j, int xi, double *Li, double (*Lijxi)[tNodes][mDomain], int (*nbr_l)[tNodes],
                        struct node *Nodes, int pLen, double *fP, double *bounds) {

    double lb_tmp, ub_tmp, eps;
    double lb, ub, x, mIter = 500, it = 0;
    double tmp[Nodes[j].domain];
    struct func_Params p;
    p.d = fP;
    param_len = pLen;


    for (int xj = 0; xj < Nodes[j].domain; ++xj) {
        tmp[xj] = -1 * Lijxi[j][i][xj];
    }
    lb_tmp = max1D(Nodes[j].domain, tmp);
    ub_tmp = Li[i] - computeSum_lijxi_k(i, j, xi, Lijxi, nbr_l);

    lb = lb_tmp;
    ub = ub_tmp;
    eps = (ub - lb) / 2;
    x = lb + eps;

    while (true) {
        if (fabs(ratFun2(x, (void *) &p)) <= zeroThreshold) {
            bounds[0] = x;
            bounds[1] = x;
            return 1;
        }
        if (ratFun2(x, (void *) &p) > 0) {
            lb = x;
            if (x != lb_tmp && ub != ub_tmp) {
                bounds[0] = x;
                bounds[1] = ub;
                return 1;
            }
            eps = (ub - lb) / 2;
            x = lb + eps;
        } else {
            ub = x;
            if (lb != lb_tmp && x != ub_tmp) {
                bounds[0] = lb;
                bounds[1] = x;
                return 1;
            }
            eps = (ub - lb) / 2;
            x = ub - eps;
        }
        it += 1;
        if (it > mIter) {
            printf("Unable to find bounds for Lijxi (%d, %d, %d) \n", i, j, xi);
            return 0;
        }
    }
}

void updateLambdaParameters(double (*PiFi)[mDomain], double *Li, double (*Lijxi)[tNodes][mDomain], struct node *Nodes,
                            double (*Pstar_ij)[tNodes][mDomain][mDomain], int lpLen, struct edge *LP,
                            int (*nbr_l)[tNodes], double (*thetaH)[tNodes][mDomain][mDomain]) {



    int liC  = 0;
    int lijxC = 0;
    int iCount = 0;
    // ----------------------- For Node i ----------------------- //
    for (int e = 0; e < lpLen; ++e) {
        int i = LP[e].i;
        int j = LP[e].j;

        double bounds[2];
        // ------------------- Update Lijxi ------------------- //
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            int pLen = 2 * Nodes[j].domain + 2;
            double fP[pLen];
            getLijxiParam(i, j, xi, fP, PiFi, Pstar_ij, Li, Lijxi, Nodes, nbr_l, thetaH);
            double bFlag = computeBounds_lijxi(i, j, xi, Li, Lijxi, nbr_l, Nodes, pLen, fP, bounds);
            double lb = bounds[0], ub = bounds[1];
            struct func_Params p;
            p.d = fP;
            param_len = pLen;

            if (lb == ub) {
                Lijxi[i][j][xi] = lb;
            } else {
                double lijxiOld = Lijxi[i][j][xi];
                if (bFlag == 1) {
                    if (!(fabs(ratFun2(lijxiOld, (void *) &p) <= zeroThreshold && lb <= lijxiOld <= ub))) {
                        double root = getRoot2(fP, pLen, lb, ub);
                        Lijxi[i][j][xi] = root;

                    }
                }
            }
        }

        // ------------------ Update Li ------------------- //
        double fP[2 * Nodes[i].domain + 1];
        int pLen = 2 * Nodes[i].domain + 1;

        bounds[0] = 0;
        bounds[1] = 0;

        for (int k = 0; k < pLen; ++k) {
            fP[k] = -1;
        }

        getLiParam(tNodes, mDomain, i, fP, PiFi, Lijxi, Nodes, nbr_l);

        getBounds_Li(tNodes, mDomain, i, pLen, fP, Lijxi, bounds, Nodes, nbr_l);
        double lb = bounds[0];
        double ub = bounds[1];
        struct func_Params p;
        p.d = fP;
        param_len = pLen;
        if (lb == ub) {
            Li[i] = lb;
        } else {
            double liOld = Li[i];
            if (!(fabs(ratFun(liOld, (void *) &p)) <= zeroThreshold && lb <= liOld <=
            ub)) {
                double root = getRoot(fP, pLen, lb, ub);
                Li[i] = root;
                // printf("New Root 2 | ");
            }
        }

        // ----------------------- For Node j ---------------------- //
        // -------------- Update Ljixj ----------------- //
        for (int xj = 0; xj < Nodes[j].domain; ++xj) {
            int pLen = 2 * Nodes[i].domain + 2;
            double fP[pLen];
            bounds[0] = 0;
            bounds[1] = 0;
            getLijxiParam(j, i, xj, fP, PiFi, Pstar_ij, Li, Lijxi, Nodes, nbr_l, thetaH);

            double bFlag = computeBounds_lijxi(j, i, xj, Li, Lijxi, nbr_l, Nodes, pLen, fP, bounds);
            double lb = bounds[0], ub = bounds[1];
            struct func_Params p;
            p.d = fP;
            param_len = pLen;

            if (lb == ub) {
                Lijxi[j][i][xj] = lb;
            } else {
                double lijxiOld = Lijxi[j][i][xj];
                if (bFlag == 1) {
                    if (!(fabs(ratFun2(lijxiOld, (void *) &p)) <= zeroThreshold && lb <= lijxiOld <= ub)) {
                        double root = getRoot2(fP, pLen, lb, ub);
                        Lijxi[j][i][xj] = root;
                    }
                }
            }
        }

        // ------------------ Update Lj ----------------- //
        double fP2[2 * Nodes[j].domain + 1];
        bounds[0] = 0;
        bounds[1] = 0;
        getLiParam(tNodes, mDomain, j, fP2, PiFi, Lijxi, Nodes, nbr_l);

        int pLen2 = 2 * Nodes[j].domain + 1;
        getBounds_Li(tNodes, mDomain, j, pLen2, fP2, Lijxi, bounds, Nodes, nbr_l);


        double lb2 = bounds[0];
        double ub2 = bounds[1];
        struct func_Params p2;
        p2.d = fP2;
        param_len = pLen2;
        if (lb2 == ub2) {
            Li[j] = lb2;

        } else {

            double liOld2 = Li[j];
            if (!(fabs(ratFun(liOld2, (void *) &p2)) <= zeroThreshold && lb2 <= liOld2 <= ub2)) {

                double root2 = getRoot(fP2, pLen2, lb2, ub2);
                Li[j] = root2;
                // printf("New Root 3 | ");

            }

        }

    }

}

void computePstar_ij(double (*Pij)[tNodes][mDomain][mDomain], double (*Pstar_ij)[tNodes][mDomain][mDomain],
                     double (*thetaH)[tNodes][mDomain][mDomain], double (*Lijxi)[tNodes][mDomain], int lpLen,
                     struct edge *LP, struct node *Nodes) {

    double Num, Den;

    for (int e = 0; e < lpLen; ++e) {
        int i = LP[e].i;
        int j = LP[e].j;
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            for (int xj = 0; xj < Nodes[j].domain; ++xj) {
                Num = thetaH[i][j][xi][xj] * Pij[i][j][xi][xj];
                Den = Lijxi[i][j][xi] + Lijxi[j][i][xj];
                Pstar_ij[i][j][xi][xj] = Num / Den;
                Pstar_ij[j][i][xj][xi] = Num / Den;
            }
        }
    }
}

double computeViolation_Pij(double (*PiFi)[mDomain], double (*Pij)[tNodes][mDomain][mDomain],
                            double (*thetaH)[tNodes][mDomain][mDomain], double *Li, double (*Lijxi)[tNodes][mDomain],
                            int lpLen, struct edge *LP, struct node *Nodes, int (*nbr_l)[tNodes]) {

    // ---------- Compute Pstar / Pstar_ij --------- //
    double Pstar[tNodes][mDomain];
    computePstar(Nodes, PiFi, Li, Lijxi, nbr_l, Pstar);

    double Pstar_ij[tNodes][tNodes][mDomain][mDomain];
    computePstar_ij(Pij, Pstar_ij, thetaH, Lijxi, lpLen, LP, Nodes);

    int nb_ij;
    double maxVio = -1e-50;
    for (int i = 0; i < tNodes; ++i) {
        for (int j = 0; j < tNodes; ++j) {
            nb_ij = nbr_l[i][j];
            if (nb_ij != -1) {
                for (int xi = 0; xi < Nodes[i].domain; ++xi) {
                    double sum = 0;
                    for (int xj = 0; xj < Nodes[j].domain; ++xj) {
                        sum += Pstar_ij[i][j][xi][xj];
                    }
                    sum -= Pstar[i][xi];

                    if (fabs(sum) > maxVio) {
                        maxVio = fabs(sum);
                    }
                }
            }
        }
    }

    if (maxVio == -1e-50) {
        return 0;
    } else {
        return maxVio;
    }
}

bool bcdConvergence(double con1, double con2){

    double rcon1 = roundDec(con1, 4);
    double rcon2 = roundDec(con2, 4);

    if(con1 <= conVioThres && con2 <= conVioThres){
        return true;
    }

    if(rcon1 == oldEst1){
        bcdConflag1 += 1;
    }
    else{
        bcdConflag1 = 0;
    }
    if(rcon2 == oldEst2){
        bcdConflag2 += 1;
    }
    else{
        bcdConflag2 = 0;
    }

    if(bcdConflag1 > bcdConv && bcdConflag2 > bcdConv){
        bcdConflag1 = 0;
        bcdConflag2 = 0;
        oldEst1 = -1;
        oldEst2 = -1;
        return true;
    }

    oldEst1 = rcon1;
    oldEst2 = rcon2;
    return false;
}

// ---------------- Newton Method ----------------- //

void updateLambda_i_newton(int n, int d, double (*PiFi)[d], double *Li, double (*Lijxi)[n][d], struct node *Nodes,
                    int (*nbr_l)[n]) {

    double root;
    double bounds[2], lb, ub;

    for (int i = 0; i < n; ++i) {
        double fP[2 * Nodes[i].domain + 1];
        bounds[0] = 0;
        bounds[1] = 0;
        getLiParam(n, d, i, fP, PiFi, Lijxi, Nodes, nbr_l);
        int pLen = 2 * Nodes[i].domain + 1;
        getBounds_Li(n, d, i, pLen, fP, Lijxi, bounds, Nodes, nbr_l);
        lb = bounds[0];
        ub = bounds[1];
        if (lb == ub) {
            Li[i] = lb;
        } else {
            // root = getRoot(fP, pLen, lb, ub);
            root = getRoot_newton(fP, pLen, lb);
            Li[i] = root;
            // printf("New Root 1 | ");
        }
    }
}

void updateLambdaParameters_newton(double (*PiFi)[mDomain], double *Li, double (*Lijxi)[tNodes][mDomain], struct node *Nodes,
                            double (*Pstar_ij)[tNodes][mDomain][mDomain], int lpLen, struct edge *LP,
                            int (*nbr_l)[tNodes], double (*thetaH)[tNodes][mDomain][mDomain]) {

    int liC  = 0;
    int lijxC = 0;
    int iCount = 0;
    // ----------------------- For Node i ----------------------- //
    for (int e = 0; e < lpLen; ++e) {
        int i = LP[e].i;
        int j = LP[e].j;

        double bounds[2];
        // ------------------- Update Lijxi ------------------- //
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            int pLen = 2 * Nodes[j].domain + 2;
            double fP[pLen];
            getLijxiParam(i, j, xi, fP, PiFi, Pstar_ij, Li, Lijxi, Nodes, nbr_l, thetaH);
            double bFlag = computeBounds_lijxi(i, j, xi, Li, Lijxi, nbr_l, Nodes, pLen, fP, bounds);
            double lb = bounds[0], ub = bounds[1];
            struct func_Params p;
            p.d = fP;
            param_len = pLen;

            if (lb == ub) {
                Lijxi[i][j][xi] = lb;
            } else {
                double lijxiOld = Lijxi[i][j][xi];
                if (bFlag == 1) {
                    if (!(fabs(ratFun2(lijxiOld, (void *) &p) <= zeroThreshold && lb <= lijxiOld <= ub))) {
                        // double root = getRoot2(fP, pLen, lb, ub);
                        double root = getRoot2_newton(fP, pLen, lb);
                        Lijxi[i][j][xi] = root;
                    }
                }
            }
        }

        // ------------------ Update Li ------------------- //
        double fP[2 * Nodes[i].domain + 1];
        int pLen = 2 * Nodes[i].domain + 1;

        bounds[0] = 0;
        bounds[1] = 0;

        for (int k = 0; k < pLen; ++k) {
            fP[k] = -1;
        }

        getLiParam(tNodes, mDomain, i, fP, PiFi, Lijxi, Nodes, nbr_l);

        getBounds_Li(tNodes, mDomain, i, pLen, fP, Lijxi, bounds, Nodes, nbr_l);
        double lb = bounds[0];
        double ub = bounds[1];
        struct func_Params p;
        p.d = fP;
        param_len = pLen;
        if (lb == ub) {
            Li[i] = lb;
        } else {
            double liOld = Li[i];
            if (!(fabs(ratFun(liOld, (void *) &p)) <= zeroThreshold && lb <= liOld <=
            ub)) {
                // double root = getRoot(fP, pLen, lb, ub);
                double root = getRoot_newton(fP, pLen, lb);
                Li[i] = root;
                // printf("New Root 2 | ");
            }
        }

        // ----------------------- For Node j ---------------------- //
        // -------------- Update Ljixj ----------------- //
        for (int xj = 0; xj < Nodes[j].domain; ++xj) {
            int pLen = 2 * Nodes[i].domain + 2;
            double fP[pLen];
            bounds[0] = 0;
            bounds[1] = 0;
            getLijxiParam(j, i, xj, fP, PiFi, Pstar_ij, Li, Lijxi, Nodes, nbr_l, thetaH);

            double bFlag = computeBounds_lijxi(j, i, xj, Li, Lijxi, nbr_l, Nodes, pLen, fP, bounds);
            double lb = bounds[0], ub = bounds[1];
            struct func_Params p;
            p.d = fP;
            param_len = pLen;

            if (lb == ub) {
                Lijxi[j][i][xj] = lb;
            } else {
                double lijxiOld = Lijxi[j][i][xj];
                if (bFlag == 1) {
                    if (!(fabs(ratFun2(lijxiOld, (void *) &p)) <= zeroThreshold && lb <= lijxiOld <= ub)) {
                        // double root = getRoot2(fP, pLen, lb, ub);
                        double root = getRoot2_newton(fP, pLen, lb);
                        Lijxi[j][i][xj] = root;
                    }
                }
            }
        }

        // ------------------ Update Lj ----------------- //
        double fP2[2 * Nodes[j].domain + 1];
        bounds[0] = 0;
        bounds[1] = 0;
        getLiParam(tNodes, mDomain, j, fP2, PiFi, Lijxi, Nodes, nbr_l);

        int pLen2 = 2 * Nodes[j].domain + 1;
        getBounds_Li(tNodes, mDomain, j, pLen2, fP2, Lijxi, bounds, Nodes, nbr_l);


        double lb2 = bounds[0];
        double ub2 = bounds[1];
        struct func_Params p2;
        p2.d = fP2;
        param_len = pLen2;
        if (lb2 == ub2) {
            Li[j] = lb2;

        } else {

            double liOld2 = Li[j];
            if (!(fabs(ratFun(liOld2, (void *) &p2)) <= zeroThreshold && lb2 <= liOld2 <= ub2)) {

                // double root2 = getRoot(fP2, pLen2, lb2, ub2);
                double root2 = getRoot_newton(fP2, pLen2, lb2);
                Li[j] = root2;
                // printf("New Root 3 | ");

            }

        }

    }

}

// ----------------------------------------------- //

void solveBCD(double (*PiFi)[mDomain], double *Li, double (*Lijxi)[tNodes][mDomain], struct node *Nodes,
              int (*nbr_l)[tNodes], double (*Pij)[tNodes][mDomain][mDomain], int lpLen, struct edge *LP,
              double (*thetaH)[tNodes][mDomain][mDomain], FILE *f) {

    double con1, con2;
    inner_iteration = 0;
    bcdRuntime = 0;
    double bcdT = 0;
    clock_t t1 = clock();

    if (iter_count == 1) {
        init_Li_Lijxi(tNodes, mDomain, nbr_l, Nodes, Lijxi, Li);
    } else {
        updateLambda_i(tNodes, mDomain, PiFi, Li, Lijxi, Nodes, nbr_l);
        updateLambdaParameters(PiFi, Li, Lijxi, Nodes, Pij, lpLen, LP, nbr_l, thetaH);
    }

    con1 = computeViolation_Pi(Nodes, PiFi, Li, Lijxi, nbr_l);
    con2 = computeViolation_Pij(PiFi, Pij, thetaH, Li, Lijxi, lpLen, LP, Nodes, nbr_l);

    clock_t t2 = clock();
    bcdT = (double)(t2 - t1) / CLOCKS_PER_SEC;
    bcdRuntime += bcdT;

    // while (roundDec(con1, 3) > conVioThres || roundDec(con2, 3) > conVioThres) {
    // while (bcdConvergence(con1, con2) == false) {

    while ((con1 > conVioThres) || (con2 > conVioThres)) {

        clock_t t1 = clock();
        inner_iteration += 1;

        // printf("\nBCD : %d ", inner_iteration);
        // fprintf(f, "\nBCD : %d ", inner_iteration);

        if (inner_iteration > 50){
            updateLambda_i_newton(tNodes, mDomain, PiFi, Li, Lijxi, Nodes, nbr_l);
            updateLambdaParameters_newton(PiFi, Li, Lijxi, Nodes, Pij, lpLen, LP, nbr_l, thetaH);
        }
        else{
        updateLambda_i(tNodes, mDomain, PiFi, Li, Lijxi, Nodes, nbr_l);
        updateLambdaParameters(PiFi, Li, Lijxi, Nodes, Pij, lpLen, LP, nbr_l, thetaH);
        }

        con1 = computeViolation_Pi(Nodes, PiFi, Li, Lijxi, nbr_l);
        con2 = computeViolation_Pij(PiFi, Pij, thetaH, Li, Lijxi, lpLen, LP, Nodes, nbr_l);

        // printf(" Con1 : %.15f", con1);
        // printf(" Con2 : %.15f", con2);


        // fprintf(f, " Con1 : %.15f", con1);
        // fprintf(f, " Con2 : %.15f", con2);


        clock_t t2 = clock();
        bcdT = (double)(t2 - t1) / CLOCKS_PER_SEC;
        bcdRuntime += bcdT;
        if (bcdRuntime > maxRuntime){
            break;
        }
    }
    printf("\n");
}

void updatePij(double (*Pij)[tNodes][mDomain][mDomain], double (*Pstar_ij)[tNodes][mDomain][mDomain]) {

    for (int i = 0; i < tNodes; ++i) {
        for (int j = 0; j < tNodes; ++j) {
            for (int xi = 0; xi < mDomain; ++xi) {
                for (int xj = 0; xj < mDomain; ++xj) {
                    Pij[i][j][xi][xj] = Pstar_ij[i][j][xi][xj];
                }
            }
        }

    }

}

void computeNbrPot(int i, int xi, double (*Pstar)[mDomain], int (*nbr_l)[tNodes], int (*nbr_q)[tNodes], struct node *Nodes,
              double (*theta)[tNodes][mDomain][mDomain], double *res) {

    double sum = 0;

    int nbr_i[tNodes];
    for (int j = 0; j < tNodes; ++j) {
        nbr_i[j] = -1;
    }

    for (int j = 0; j < tNodes; ++j) {
        if (nbr_l[i][j] != -1) {
            nbr_i[j] = nbr_l[i][j];
        }
        if (nbr_q[i][j] != -1) {
            nbr_i[j] = nbr_q[i][j];
        }
    }

    for (int j = 0; j < tNodes; ++j) {
        if (nbr_i[j] != -1) {
            for (int xj = 0; xj < Nodes[j].domain; ++xj) {
                sum += Pstar[j][xj] * theta[i][j][xi][xj];
            }
        }
    }

    res[0] = (double) xi;
    res[1] = sum;

}

int compareArrays(double *a, double *b, int n) {

    int ii;


    for (ii = 0; ii < n; ++ii) {
        printf("%f %f\n", a[ii], b[ii]);
    }

    printf("----------\n");


    for (ii = 0; ii < n; ++ii) {
        if (a[ii] != b[ii]) {
            return 0;
        }
    }
    return 1;
}

bool primalConvergence(int (*tmpSol)[tNodes], int it_Count) {

    int cFlag1, cFlag2;

    if (it_Count >= 2) {

        for (int k = 0; k < tNodes; ++k) {
            if (tmpSol[0][k] != tmpSol[1][k]) {
                return true;
            }
        }

        for (int k = 0; k < tNodes; ++k) {
            if (tmpSol[1][k] != tmpSol[2][k]) {
                return true;
            }
        }

        return false;
    } else {
        return true;
    }

}

void primalExtract(double (*Pstar)[mDomain], int *sol, int (*nbr_l)[tNodes], int (*nbr_q)[tNodes], struct node *Nodes,
                   double (*theta)[tNodes][mDomain][mDomain]) {

    int it_Count = 0;
    int tmpSol[3][tNodes], tSolCount = 0;
    double Pstar_tmp[tNodes][mDomain];

    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < tNodes; ++i) {
            tmpSol[j][i] = -1;
        }
    }

    for (int l = 0; l < tNodes; ++l) {
        for (int xt = 0; xt < Nodes[l].domain; ++xt) {
            Pstar_tmp[l][xt] = Pstar[l][xt];
        }
    }

    bool conFlag = true;

    while (conFlag) {
        for (int i = 0; i < tNodes; ++i) {
            double maxm[2] = {-1, -1e-100};
            for (int xi = 0; xi < Nodes[i].domain; ++xi) {
                double temp[2];
                computeNbrPot(i, xi, Pstar_tmp, nbr_l, nbr_q, Nodes, theta, temp);
                if (maxm[1] < temp[1]) {
                    maxm[0] = temp[0];
                    maxm[1] = temp[1];
                }
            }
            sol[i] = (int) maxm[0];

            // ----- Update Pstar Temp ------ //
            for (int xtmp = 0; xtmp < Nodes[i].domain; ++xtmp) {
                if (xtmp == sol[i]) {
                    Pstar_tmp[i][xtmp] = 1.0;
                } else {
                    Pstar_tmp[i][xtmp] = 0.0;
                }
            }
        }

        double one = 1;
        for (int i = 0; i < tNodes; ++i) {
            for (int xi = 0; xi < Nodes[i].domain; ++xi) {
                if (Pstar_tmp[i][xi] == one)
                {
                    tmpSol[tSolCount][i] = xi;
                }
            }
        }

        conFlag = primalConvergence(tmpSol, it_Count);

        it_Count += 1;
        tSolCount += 1;
        if (tSolCount >= 3) {
            tSolCount = 0;
        }

        if (it_Count > 500) {
            printf("Primal Extract Not Converged \n");
            printf("Exiting !\n");
            halt();
        }

    }

}

void checkThreshold_P(struct node *Nodes, double (*Pstar)[mDomain]) {

    for (int i = 0; i < tNodes; ++i) {
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            if (Pstar[i][xi] < P_Threshold) {
                Pstar[i][xi] = P_Threshold;
            }
        }
    }
}

void checkThreshold_Pij(double (*Pij)[tNodes][mDomain][mDomain]) {

    for (int i = 0; i < tNodes; ++i) {
        for (int j = 0; j < tNodes; ++j) {
            for (int xi = 0; xi < mDomain; ++xi) {
                for (int xj = 0; xj < mDomain; ++xj) {
                    if (Pij[i][j][xi][xj] < P_Threshold) {
                        Pij[i][j][xi][xj] = P_Threshold;
                    }

                }
            }
        }
    }

}

double computeIntObjective(double (*pot)[tNodes][mDomain][mDomain], int qpLen, struct edge *QP, int lpLen, struct edge *LP, int *sol){

    double tmpQP = 0, tmpLP = 0;
    for (int e = 0; e < qpLen; ++e) {
        int i = QP[e].i;
        int j = QP[e].j;
        int xi = sol[i];
        int xj = sol[j];
        tmpQP += pot[i][j][xi][xj];
    }

    for (int e = 0; e < lpLen; ++e) {
        int i = LP[e].i;
        int j = LP[e].j;
        int xi = sol[i];
        int xj = sol[j];
        tmpLP += pot[i][j][xi][xj];
    }

    return tmpQP + tmpLP;
}

double computeObjective(double (*pot)[tNodes][mDomain][mDomain], double (*Pstar)[mDomain], double (*Pstar_ij)[tNodes][mDomain][mDomain], int qpLen, struct edge *QP, int lpLen, struct edge *LP, struct node *Nodes) {

    double tmpQP = 0, tmpLP = 0;
    for (int e = 0; e < qpLen; ++e) {
        int i = QP[e].i;
        int j = QP[e].j;
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            for (int xj = 0; xj < Nodes[j].domain; ++xj) {
                double Pixi = Pstar[i][xi];
                double Pjxj = Pstar[j][xj];
                if (Pixi <= P_Threshold)
                    Pixi = 0.0;
                if (Pjxj <= P_Threshold)
                    Pjxj = 0.0;
                tmpQP += Pixi * Pjxj * pot[i][j][xi][xj];
            }
        }
    }

    for (int e = 0; e < lpLen; ++e) {
        int i = LP[e].i;
        int j = LP[e].j;
        for (int xi = 0; xi < Nodes[i].domain; ++xi) {
            for (int xj = 0; xj < Nodes[j].domain; ++xj) {
                double Pijxixj = Pstar_ij[i][j][xi][xj];
                if (Pijxixj <= P_Threshold) {
                    Pijxixj = 0.0;
                }
                tmpLP += Pijxixj * pot[i][j][xi][xj];
            }
        }
    }

    return tmpQP + tmpLP;
}

void test2(double (*PiFi)[mDomain]) {
    printf("%f \n", PiFi[0][0]);
    printf("%d \n", mDomain);
    printf("%d \n", tNodes);

}

#endif //CLIONPROJECTS_RCDCOP_H

