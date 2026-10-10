#include "RTK_Structs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Eigen/Dense>

/*------------------------------------------------------------------------------
* lambda.cpp : integer ambiguity resolution (LAMBDA / MLAMBDA method)
*
* Reference:
*   [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*       a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82, 1995
*   [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*       integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*-----------------------------------------------------------------------------*/

#define LOOPMAX     500000          /* 最大搜索循环计数 */
#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

// 内部轻量级动态矩阵辅助函数
static double *mat(int n, int m)
{
    if (n <= 0 || m <= 0) return NULL;
    return (double *)malloc(sizeof(double) * n * m);
}

static double *zeros(int n, int m)
{
    if (n <= 0 || m <= 0) return NULL;
    return (double *)calloc(n * m, sizeof(double));
}

static double *eye(int n)
{
    double *p = zeros(n, n);
    if (p) {
        for (int i = 0; i < n; i++) p[i + i * n] = 1.0;
    }
    return p;
}

// 矩阵线性方程组求解辅助: 求解 Z^T * F = E (利用工程已有的 Eigen 求解，无需外部 BLAS/LAPACK 库)
static int solve_mat(const char *tr, const double *A, const double *Y, int n, int m, double *X)
{
    Eigen::Map<const Eigen::MatrixXd> matA(A, n, n);
    Eigen::Map<const Eigen::MatrixXd> matY(Y, n, m);
    Eigen::Map<Eigen::MatrixXd> matX(X, n, m);
    if (tr[0] == 'T' || tr[0] == 't') {
        matX = matA.transpose().partialPivLu().solve(matY);
    } else {
        matX = matA.partialPivLu().solve(matY);
    }
    return 0;
}

/* LD 分解 (Q = L' * diag(D) * L) -------------------------------------------*/
int LD(int n, const double *Q, double *L, double *D)
{
    int i, j, k, info = 0;
    double a, *A = mat(n, n);
    if (!A) return -1;
    
    memcpy(A, Q, sizeof(double) * n * n);
    for (i = n - 1; i >= 0; i--) {
        if ((D[i] = A[i + i * n]) <= 0.0) { info = -1; break; }
        a = sqrt(D[i]);
        for (j = 0; j <= i; j++) L[i + j * n] = A[i + j * n] / a;
        for (j = 0; j <= i - 1; j++) {
            for (k = 0; k <= j; k++) A[j + k * n] -= L[i + k * n] * L[i + j * n];
        }
        for (j = 0; j <= i; j++) L[i + j * n] /= L[i + i * n];
    }
    free(A);
    if (info) fprintf(stderr, "%s : LD factorization error\n", __FILE__);
    return info;
}

/* 整数高斯变换 --------------------------------------------------------------*/
void gauss(int n, double *L, double *Z, int i, int j)
{
    int k, mu;
    if ((mu = (int)ROUND(L[i + j * n])) != 0) {
        for (k = i; k < n; k++) L[k + n * j] -= (double)mu * L[k + i * n];
        for (k = 0; k < n; k++) Z[k + n * j] -= (double)mu * Z[k + i * n];
    }
}

/* 条件方差排序与列置换 ------------------------------------------------------*/
void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta, lam, a0, a1;
    
    eta = D[j] / del;
    lam = D[j + 1] * L[j + 1 + j * n] / del;
    D[j] = eta * D[j + 1];
    D[j + 1] = del;
    for (k = 0; k <= j - 1; k++) {
        a0 = L[j + k * n];
        a1 = L[j + 1 + k * n];
        L[j + k * n]     = -L[j + 1 + j * n] * a0 + a1;
        L[j + 1 + k * n] = eta * a0 + lam * a1;
    }
    L[j + 1 + j * n] = lam;
    for (k = j + 2; k < n; k++) SWAP(L[k + j * n], L[k + (j + 1) * n]);
    for (k = 0; k < n; k++)     SWAP(Z[k + j * n], Z[k + (j + 1) * n]);
}

/* LAMBDA 降相关（整数 Z 变换） ----------------------------------------------*/
void reduction(int n, double *L, double *D, double *Z)
{
    int i, j, k;
    double del;
    
    j = n - 2;
    k = n - 2;
    while (j >= 0) {
        if (j <= k) {
            for (i = j + 1; i < n; i++) gauss(n, L, Z, i, j);
        }
        del = D[j] + L[j + 1 + j * n] * L[j + 1 + j * n] * D[j + 1];
        if (del + 1E-6 < D[j + 1]) { // 考虑浮点数值误差
            perm(n, L, D, j, del, Z);
            k = j;
            j = n - 2;
        } else {
            j--;
        }
    }
}

/* 修正的收缩椭球搜索 (MLAMBDA 算法) -----------------------------------------*/
int search(int n, int m, const double *L, const double *D,
           const double *zs, double *zn, double *s)
{
    int i, j, k, c, nn = 0, imax = 0;
    double newdist, maxdist = 1E99, y;
    double *S = zeros(n, n);
    double *dist = mat(n, 1);
    double *zb = mat(n, 1);
    double *z = mat(n, 1);
    double *step = mat(n, 1);
    if (!S || !dist || !zb || !z || !step) {
        free(S); free(dist); free(zb); free(z); free(step);
        return -1;
    }
    
    k = n - 1;
    dist[k] = 0.0;
    zb[k] = zs[k];
    z[k] = ROUND(zb[k]);
    y = zb[k] - z[k];
    step[k] = SGN(y);

    for (c = 0; c < LOOPMAX; c++) {
        newdist = dist[k] + y * y / D[k];
        if (newdist < maxdist) {
            if (k != 0) {
                dist[--k] = newdist;
                for (i = 0; i <= k; i++) {
                    S[k + i * n] = S[k + 1 + i * n] + (z[k + 1] - zb[k + 1]) * L[k + 1 + i * n];
                }
                zb[k] = zs[k] + S[k + k * n];
                z[k] = ROUND(zb[k]);
                y = zb[k] - z[k];
                step[k] = SGN(y);
            } else {
                if (nn < m) {
                    if (nn == 0 || newdist > s[imax]) imax = nn;
                    for (i = 0; i < n; i++) zn[i + nn * n] = z[i];
                    s[nn++] = newdist;
                } else {
                    if (newdist < s[imax]) {
                        for (i = 0; i < n; i++) zn[i + imax * n] = z[i];
                        s[imax] = newdist;
                        for (i = imax = 0; i < m; i++) if (s[imax] < s[i]) imax = i;
                    }
                    maxdist = s[imax];
                }
                z[0] += step[0];
                y = zb[0] - z[0];
                step[0] = -step[0] - SGN(step[0]);
            }
        } else {
            if (k == n - 1) {
                break;
            } else {
                k++;
                z[k] += step[k];
                y = zb[k] - z[k];
                step[k] = -step[k] - SGN(step[k]);
            }
        }
    }

    // 按残差平方和升序排序 (s[0] <= s[1])
    for (i = 0; i < m - 1; i++) {
        for (j = i + 1; j < m; j++) {
            if (s[i] < s[j]) continue;
            SWAP(s[i], s[j]);
            for (k = 0; k < n; k++) SWAP(zn[k + i * n], zn[k + j * n]);
        }
    }

    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c >= LOOPMAX) {
        fprintf(stderr, "%s : search loop count overflow\n", __FILE__);
        return -1;
    }
    return 0;
}

/* LAMBDA 顶层接口 -----------------------------------------------------------
* 参数说明：
*   n   : 待估双差模糊度维数
*   m   : 输出候选组数 (通常为 2，即最优与次优解)
*   a   : 浮点模糊度向量 (n x 1)
*   Q   : 浮点模糊度协方差矩阵 (n x n，列优先存储)
*   F   : 输出固定整数模糊度 (n x m)
*   s   : 输出各候选解残差平方和 (1 x m)
* 返回值：0 表示成功，非 0 表示失败
*-----------------------------------------------------------------------------*/
int lambda(int n, int m, const double *a, const double *Q, double *F, double *s)
{
    int info = -1;
    double *L, *D, *Z, *z, *E;
    
    if (n <= 0 || m <= 0) return -1;
    L = zeros(n, n);
    D = mat(n, 1);
    Z = eye(n);
    z = mat(n, 1);
    E = mat(n, m);
    if (!L || !D || !Z || !z || !E) {
        free(L); free(D); free(Z); free(z); free(E);
        return -1;
    }
    
    // 1. LD 分解: Q = L^T * D * L
    if (!(info = LD(n, Q, L, D))) {
        // 2. 降相关处理: 获得整数变换矩阵 Z 及对角化更好的 L 与 D
        reduction(n, L, D, Z);
        
        // 3. 模糊度浮点解正变换: z = Z^T * a
        {
            Eigen::Map<const Eigen::MatrixXd> matZ(Z, n, n);
            Eigen::Map<const Eigen::VectorXd> veca(a, n);
            Eigen::Map<Eigen::VectorXd> vecz(z, n);
            vecz = matZ.transpose() * veca;
        }
        
        // 4. MLAMBDA 椭球搜索: 搜索 z 空间的最优 m 组整数解 E
        if (!(info = search(n, m, L, D, z, E, s))) {
            // 5. 整数解逆变换回原模糊度空间: 求解 Z^T * F = E 即 F = Z^{-T} * E
            info = solve_mat("T", Z, E, n, m, F);
        }
    }
    
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
