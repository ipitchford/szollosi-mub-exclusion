/*
 * kernel.c -- trusted interval kernel of the INDEPENDENT CHECKER (no code shared with the generator).
 *
 * Arithmetic model.  Every floating-point operation is IEEE-754 binary64 round-to-nearest (the default mode;
 * the code never changes the rounding mode), and every computed endpoint is then moved ONE ulp outward
 * (up1 / dn1).  Since a correctly rounded result lies within half a spacing of the exact result on either side,
 * and the neighbouring double lies at least that far away, dn1(fl(x op y)) <= x op y <= up1(fl(x op y)) for
 * x, y doubles and op in {+,-,*,/,sqrt}.  Interval operations are the textbook inf-sup operations built from
 * these.  No library function other than sqrt is used on any enclosure; there is no BLAS and no trigonometry.
 * Correctness therefore rests only on IEEE correct rounding of + - * / sqrt (and on the table, verified
 * separately with mpmath interval arithmetic before it is passed in).
 *
 * Angles.  All angles are integers on the grid G17 = (2 pi / 2^17) Z ("units").  The table T[m] encloses
 * exp(2 pi i m / 2^17) with componentwise error <= 2^-53.  A = 2 pi / 2^17 = pi / 2^16.
 *
 * Enclosure of exp(i phi) for phi in [c - w, c + w] (units):  exp(i c) * exp(i t), |t| <= w A, with
 *   Re exp(it) in [cos(wA), 1]  (valid for wA <= pi, i.e. w <= 2^16; else [-1, 1])
 *   Im exp(it) in [-sin(wA), sin(wA)]  (valid for wA <= pi/2, i.e. w <= 2^15; else [-1, 1])
 * where cos(wA), sin(wA) are enclosed by the table entry T[w].
 */
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define NG 131072
#define NGMASK (NG - 1)
#define MAXP 4
#define MAXT 8
#define MAXN 7

typedef struct { double lo, hi; } iv;
typedef struct { iv re, im; } civ;

/* ---------------------------------------------------------------- one-ulp outward steps */
static inline double up1(double x) {
    if (x != x) return x;                      /* NaN propagates */
    if (x == INFINITY) return x;
    if (x == 0.0) return 4.9406564584124654e-324;
    union { double d; int64_t i; } u; u.d = x;
    if (x > 0) u.i += 1; else u.i -= 1;
    return u.d;
}
static inline double dn1(double x) { return -up1(-x); }

/* ---------------------------------------------------------------- real intervals */
static const iv IV_ALL = { -INFINITY, INFINITY };
static inline iv ivp(double a) { iv r = { a, a }; return r; }
#define IVZ(a) ((a).lo == 0.0 && (a).hi == 0.0)      /* the exact zero interval (its term contributes nothing) */
static inline iv iv_add(iv a, iv b) { iv r = { dn1(a.lo + b.lo), up1(a.hi + b.hi) }; return r; }
static inline iv iv_sub(iv a, iv b) { iv r = { dn1(a.lo - b.hi), up1(a.hi - b.lo) }; return r; }
static inline iv iv_neg(iv a) { iv r = { -a.hi, -a.lo }; return r; }
static inline iv iv_mul(iv a, iv b) {
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    if (p1 != p1 || p2 != p2 || p3 != p3 || p4 != p4) return IV_ALL;
    double mn = p1, mx = p1;
    if (p2 < mn) mn = p2; if (p2 > mx) mx = p2;
    if (p3 < mn) mn = p3; if (p3 > mx) mx = p3;
    if (p4 < mn) mn = p4; if (p4 > mx) mx = p4;
    iv r = { dn1(mn), up1(mx) }; return r;
}
static inline iv iv_scal(double s, iv a) {         /* exact double s times interval */
    double p1 = s * a.lo, p2 = s * a.hi;
    if (p1 != p1 || p2 != p2) return IV_ALL;
    iv r; if (p1 <= p2) { r.lo = dn1(p1); r.hi = up1(p2); } else { r.lo = dn1(p2); r.hi = up1(p1); }
    return r;
}
static inline iv iv_sqr(iv a) {
    double l2 = a.lo * a.lo, h2 = a.hi * a.hi;
    iv r;
    if (a.lo >= 0) { r.lo = dn1(l2); r.hi = up1(h2); }
    else if (a.hi <= 0) { r.lo = dn1(h2); r.hi = up1(l2); }
    else { r.lo = 0.0; r.hi = up1(l2 > h2 ? l2 : h2); }
    if (r.lo < 0) r.lo = 0.0;
    if (r.hi != r.hi) return IV_ALL;
    return r;
}
static inline double iv_mag(iv a) { double x = fabs(a.lo), y = fabs(a.hi); return x > y ? x : y; }
static inline iv iv_sym_mul(iv a, double R) {      /* a * [-R, R], R >= 0 */
    double m = up1(iv_mag(a) * R);
    if (m != m) m = INFINITY;
    iv r = { -m, m }; return r;
}

/* ---------------------------------------------------------------- complex rectangles */
static inline civ c_add(civ a, civ b) { civ r = { iv_add(a.re, b.re), iv_add(a.im, b.im) }; return r; }
static inline civ c_mul(civ a, civ b) {
    civ r = { iv_sub(iv_mul(a.re, b.re), iv_mul(a.im, b.im)), iv_add(iv_mul(a.re, b.im), iv_mul(a.im, b.re)) };
    return r;
}
static inline civ c_muli(civ a) { civ r = { iv_neg(a.im), a.re }; return r; }       /* i a (exact) */
static inline civ c_neg(civ a) { civ r = { iv_neg(a.re), iv_neg(a.im) }; return r; }
static inline civ c_scal(double s, civ a) { civ r = { iv_scal(s, a.re), iv_scal(s, a.im) }; return r; }
static inline iv c_re_conj_mul(civ a, civ b) {      /* Re(conj(a) b) = a.re b.re + a.im b.im */
    return iv_add(iv_mul(a.re, b.re), iv_mul(a.im, b.im));
}
static inline iv c_abs2(civ a) { return iv_add(iv_sqr(a.re), iv_sqr(a.im)); }
static const civ C_ZERO = { { 0.0, 0.0 }, { 0.0, 0.0 } };
/* acc + e * a for an integer e: e = 0 contributes nothing (exactly), e = +-1 is an exact sign change */
static inline civ c_axpy_int(civ acc, int e, civ a) {
    if (e == 0) return acc;
    if (e == 1) return c_add(acc, a);
    if (e == -1) return c_add(acc, c_neg(a));
    return c_add(acc, c_scal((double)e, a));
}

/* ---------------------------------------------------------------- table and family (globals, set once) */
static civ *TR = 0;                 /* TR[m] rectangle enclosing exp(2 pi i m / NG) */
static double A_UP;                 /* >= pi / 2^16 */
static double TR_EPS = 1.1102230246251565e-16;   /* 2^-53 */

static int FM;                      /* number of parameters */
static int FE[6][6][MAXP];          /* exponents E_kj */
static civ FCC[6][6];               /* conj(coef_kj) rectangles */
static int FCK[6][6];               /* 1: coef == +1 exactly, -1: coef == -1 exactly, 0: general */
static int FNT;                     /* constraint terms (0: no constraint) */
static int FCE[MAXT][MAXP];
static civ FCT[MAXT];               /* constraint coefficients */

int set_table(const double *T, double pi_up) {
    if (TR) free(TR);
    TR = (civ *)malloc(sizeof(civ) * NG);
    if (!TR) return -1;
    for (int m = 0; m < NG; m++) {
        double re = T[2 * m], im = T[2 * m + 1];
        TR[m].re.lo = dn1(re - TR_EPS); TR[m].re.hi = up1(re + TR_EPS);
        TR[m].im.lo = dn1(im - TR_EPS); TR[m].im.hi = up1(im + TR_EPS);
    }
    A_UP = pi_up / 65536.0;          /* exact scaling by a power of two */
    return 0;
}

int set_family(int m, const int32_t *E, const double *coef_rect, int nt, const int32_t *CE, const double *ct_rect) {
    if (m > MAXP || nt > MAXT) return -1;
    FM = m; FNT = nt;
    for (int k = 0; k < 6; k++) for (int j = 0; j < 6; j++) {
        for (int i = 0; i < m; i++) FE[k][j][i] = E[(k * 6 + j) * m + i];
        const double *c = coef_rect + 4 * (k * 6 + j);      /* coef rectangle: re_lo re_hi im_lo im_hi */
        civ cc; cc.re.lo = c[0]; cc.re.hi = c[1]; cc.im.lo = -c[3]; cc.im.hi = -c[2];   /* conjugate */
        FCC[k][j] = cc;
        if (c[0] == 1.0 && c[1] == 1.0 && c[2] == 0.0 && c[3] == 0.0) FCK[k][j] = 1;
        else if (c[0] == -1.0 && c[1] == -1.0 && c[2] == 0.0 && c[3] == 0.0) FCK[k][j] = -1;
        else FCK[k][j] = 0;
    }
    for (int t = 0; t < nt; t++) {
        for (int i = 0; i < m; i++) FCE[t][i] = CE[t * m + i];
        const double *c = ct_rect + 4 * t;
        FCT[t].re.lo = c[0]; FCT[t].re.hi = c[1]; FCT[t].im.lo = c[2]; FCT[t].im.hi = c[3];
    }
    return 0;
}

static inline civ tab(int64_t c) { return TR[c & NGMASK]; }         /* two's complement: c mod NG */

/* rectangle enclosing exp(i t), |t| <= w A, w >= 0 units */
static inline civ rot(int64_t w) {
    civ r;
    if (w <= 0) { r.re.lo = 1.0; r.re.hi = 1.0; r.im.lo = 0.0; r.im.hi = 0.0; return r; }
    double cl = (w <= NG / 2) ? TR[w].re.lo : -1.0;
    double sm = (w <= NG / 4) ? TR[w].im.hi : 1.0;
    if (cl < -1.0) cl = -1.0;
    if (sm > 1.0) sm = 1.0;
    r.re.lo = cl; r.re.hi = 1.0; r.im.lo = -sm; r.im.hi = sm;
    return r;
}
/* a * [cl, 1] for 0 <= cl <= 1 (exact endpoints where the factor is 1) */
static inline iv iv_mul_cl1(iv a, double cl) {
    iv r;
    if (a.lo >= 0) { r.lo = dn1(a.lo * cl); r.hi = a.hi; }
    else if (a.hi <= 0) { r.lo = a.lo; r.hi = up1(a.hi * cl); }
    else { r.lo = a.lo; r.hi = a.hi; }
    return r;
}
/* enclosure of exp(i phi), phi in [c - w, c + w] units.
 * Fast path (w <= NG/4, so rot(w) = [cl, 1] x [-s, s] with cl >= 0): the rectangle product
 *   (A + iB)(R + iS), R in [cl, 1], S in [-s, s]:  Re = A R - B S,  Im = B R + A S,
 * where A R and B R use iv_mul_cl1 and B S, A S are symmetric intervals [-|B| s, |B| s].  Identical in value to
 * the general c_mul (the same inf-sup products), only cheaper. */
static inline civ expi(int64_t c, int64_t w) {
    if (w <= 0) return tab(c);
    if (w <= NG / 4) {
        civ a = tab(c);
        double cl = TR[w].re.lo, sm = TR[w].im.hi;
        if (cl < 0.0) return c_mul(a, rot(w));
        if (sm > 1.0) sm = 1.0;
        civ r;
        r.re = iv_sub(iv_mul_cl1(a.re, cl), iv_sym_mul(a.im, sm));
        r.im = iv_add(iv_mul_cl1(a.im, cl), iv_sym_mul(a.re, sm));
        return r;
    }
    return c_mul(tab(c), rot(w));
}
static inline double angle_up(int64_t n) {           /* >= n * pi / 2^16 */
    return up1((double)n * A_UP);
}

/* ---------------------------------------------------------------- terms of S_j = sum_k conj(H_kj) v_k
 * term(k, j) = conj(coef_kj) exp(i (theta_k - <E_kj, p>)), theta_1 = 0.
 * c[k], n[k] (k = 0..4): centre and half-width (units) of theta_{k+2}.
 * pc[i]: parameter centres (units); ph: parameter half-width (units); vt: theta varies; vp[i]: p_i varies. */
static inline civ term(int k, int j, const int64_t *c, const int64_t *n, const int64_t *pc, int64_t ph,
                       int vt, const int *vp) {
    int64_t cc = (k == 0) ? 0 : c[k - 1];
    int64_t w = (k == 0 || !vt) ? 0 : n[k - 1];
    for (int i = 0; i < FM; i++) {
        int e = FE[k][j][i];
        cc -= (int64_t)e * pc[i];
        if (vp[i]) w += (int64_t)(e < 0 ? -e : e) * ph;
    }
    civ t = expi(cc, w);
    int kind = FCK[k][j];
    if (kind == 1) return t;
    if (kind == -1) return c_neg(t);
    return c_mul(FCC[k][j], t);
}
static inline civ cterm(int t_, const int64_t *pc, int64_t ph, const int *vp) {  /* constraint term */
    int64_t cc = 0, w = 0;
    for (int i = 0; i < FM; i++) {
        int e = FCE[t_][i];
        cc += (int64_t)e * pc[i];
        if (vp[i]) w += (int64_t)(e < 0 ? -e : e) * ph;
    }
    return c_mul(FCT[t_], expi(cc, w));
}

/* ---------------------------------------------------------------- modulus test
 * Cheap disc form: S_j lies within dist R_j of the rectangle M_j = sum_k conj(coef) T[c_kj] (all errors of
 * the table and coefficients inside M_j), R_j = sum_k w_kj A (chord <= arc).  Excluded iff for some j
 * max|S_j| < sqrt6  or  min|S_j| > sqrt6  (compared as squares). */
static int modulus_core(int64_t nb, const int64_t *lo, const int64_t *hi, const int64_t *idx, int64_t s,
                        const int64_t *pc, int64_t ph, uint8_t *excl) {
    /* parameter contributions do not depend on the phase box: cc_p = -<E_kj, pc>, w_p = sum |E_kj,i| ph */
    int64_t ccp[6][6], wp[6][6];
    for (int k = 0; k < 6; k++) for (int j = 0; j < 6; j++) {
        int64_t cc = 0, w = 0;
        for (int i = 0; i < FM; i++) { int e = FE[k][j][i]; cc -= (int64_t)e * pc[i]; w += (int64_t)(e < 0 ? -e : e) * ph; }
        ccp[k][j] = cc; wp[k][j] = w;
    }
    for (int64_t b = 0; b < nb; b++) {
        int64_t c[5], n[5];
        if (idx) { for (int k = 0; k < 5; k++) { c[k] = 2 * idx[5 * b + k] * s + s; n[k] = s; } }
        else { for (int k = 0; k < 5; k++) { c[k] = lo[5 * b + k] + hi[5 * b + k]; n[k] = hi[5 * b + k] - lo[5 * b + k]; } }
        uint8_t ex = 0;
        for (int j = 0; j < 6 && !ex; j++) {
            civ M = C_ZERO; int64_t W = 0;
            for (int k = 0; k < 6; k++) {
                int64_t cc = ccp[k][j] + ((k == 0) ? 0 : c[k - 1]);
                W += wp[k][j] + ((k == 0) ? 0 : n[k - 1]);
                civ t = tab(cc);
                int kind = FCK[k][j];
                if (kind == -1) t = c_neg(t); else if (kind == 0) t = c_mul(FCC[k][j], t);
                M = c_add(M, t);
            }
            double R = angle_up(W);                    /* sum of chord bounds */
            /* |S| <= |M|max + R,  |S| >= |M|min - R */
            iv mre = M.re, mim = M.im;
            double maxre = iv_mag(mre), maxim = iv_mag(mim);
            double minre = (mre.lo > 0) ? mre.lo : ((mre.hi < 0) ? -mre.hi : 0.0);
            double minim = (mim.lo > 0) ? mim.lo : ((mim.hi < 0) ? -mim.hi : 0.0);
            double Mmax = up1(sqrt(up1(up1(maxre * maxre) + up1(maxim * maxim))));
            double Mmin = dn1(sqrt(dn1(dn1(minre * minre) + dn1(minim * minim))));
            double smax = up1(Mmax + R);
            double smin = dn1(Mmin - R);
            if (up1(smax * smax) < 6.0) ex = 1;
            else if (smin > 0 && dn1(smin * smin) > 6.0) ex = 1;
        }
        excl[b] = ex;
    }
    return 0;
}
int modulus_test(int64_t nb, const int64_t *lo, const int64_t *hi, const int64_t *pc, int64_t ph, uint8_t *excl) {
    return modulus_core(nb, lo, hi, 0, 0, pc, ph, excl);
}
/* same test for dyadic boxes given by index and level L (G16 bounds lo = idx 2^(16-L), hi = lo + 2^(16-L)) */
int modulus_test_idx(int64_t nb, const int64_t *idx, int L, const int64_t *pc, int64_t ph, uint8_t *excl) {
    return modulus_core(nb, 0, 0, idx, (int64_t)1 << (16 - L), pc, ph, excl);
}

/* ---------------------------------------------------------------- dense small linear algebra (numerics only) */
static int invert(int n, double a[MAXN][MAXN], double y[MAXN][MAXN]) {
    double m[MAXN][2 * MAXN];
    for (int i = 0; i < n; i++) for (int j = 0; j < 2 * n; j++) m[i][j] = (j < n) ? a[i][j] : (j - n == i ? 1.0 : 0.0);
    for (int col = 0; col < n; col++) {
        int piv = col; double best = fabs(m[col][col]);
        for (int r = col + 1; r < n; r++) if (fabs(m[r][col]) > best) { best = fabs(m[r][col]); piv = r; }
        if (!(best > 1e-300)) return -1;
        if (piv != col) for (int j = 0; j < 2 * n; j++) { double t = m[col][j]; m[col][j] = m[piv][j]; m[piv][j] = t; }
        double d = m[col][col];
        for (int j = 0; j < 2 * n; j++) m[col][j] /= d;
        for (int r = 0; r < n; r++) if (r != col) {
            double f = m[r][col];
            if (f != 0.0) for (int j = 0; j < 2 * n; j++) m[r][j] -= f * m[col][j];
        }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { y[i][j] = m[i][n + j]; if (!isfinite(y[i][j])) return -1; }
    return 0;
}
static inline double midp(iv a) { return 0.5 * a.lo + 0.5 * a.hi; }

/* ---------------------------------------------------------------- parametric Krawczyk exclusion test
 * Unknowns w = (theta_2..theta_6, p_dep...) (n = 5 + nd), parameters q = remaining angles (nq = m - nd).
 * Equations: F_j = |S_j|^2 - 6 (j = 2..6) and, if the family has a constraint, Re C, Im C.
 * K - w~ = -Y [F(w~; c_q) + Fq(w~; Q)(q - c_q)] + (I - Y [J(W; c_q) + Jq(W; Q)(q - c_q)]) (W - w~)
 * (mean-value forms in q; any real Y).  Excluded iff some coordinate of K - w~ misses [-r_i, r_i]. */
static int krawczyk_core(int64_t nb, const int64_t *lo, const int64_t *hi, const int64_t *pc, int64_t ph,
                  int nd, const int32_t *dep, uint8_t *excl, double *kd_out, double *y_out) {
    int n = 5 + nd;
    int free_[MAXP], nq = 0, isdep[MAXP] = {0, 0, 0, 0};
    for (int i = 0; i < nd; i++) isdep[dep[i]] = 1;
    for (int i = 0; i < FM; i++) if (!isdep[i]) free_[nq++] = i;
    int v_none[MAXP] = {0, 0, 0, 0}, v_q[MAXP] = {0, 0, 0, 0}, v_d[MAXP] = {0, 0, 0, 0}, v_all[MAXP] = {1, 1, 1, 1};
    for (int i = 0; i < FM; i++) { v_q[i] = !isdep[i]; v_d[i] = isdep[i]; }
    double hq = angle_up(ph);                              /* parameter half-width (radians, upper bound) */
    /* constraint quantities do not depend on the phase box: compute once */
    civ C0 = C_ZERO, dCq[MAXP], dCd[MAXP], d2Cqd[MAXP][MAXP];
    for (int l = 0; l < MAXP; l++) { dCq[l] = C_ZERO; dCd[l] = C_ZERO; for (int i = 0; i < MAXP; i++) d2Cqd[l][i] = C_ZERO; }
    double dCd_mid[MAXP][2];
    if (nd) {
        for (int t = 0; t < FNT; t++) {
            C0 = c_add(C0, cterm(t, pc, ph, v_none));
            civ tq = cterm(t, pc, ph, v_q), td = cterm(t, pc, ph, v_d), ta = cterm(t, pc, ph, v_all);
            for (int l = 0; l < nq; l++) dCq[l] = c_axpy_int(dCq[l], FCE[t][free_[l]], c_muli(tq));
            for (int i = 0; i < nd; i++) dCd[i] = c_axpy_int(dCd[i], FCE[t][dep[i]], c_muli(td));
            for (int l = 0; l < nq; l++) for (int i = 0; i < nd; i++)
                d2Cqd[l][i] = c_axpy_int(d2Cqd[l][i], -FCE[t][free_[l]] * FCE[t][dep[i]], ta);
        }
        for (int i = 0; i < nd; i++) {
            civ z = C_ZERO;
            for (int t = 0; t < FNT; t++) z = c_axpy_int(z, FCE[t][dep[i]], c_muli(cterm(t, pc, ph, v_none)));
            dCd_mid[i][0] = midp(z.re); dCd_mid[i][1] = midp(z.im);
        }
    }
    for (int64_t b = 0; b < nb; b++) {
        int64_t c[5], nn[5];
        for (int k = 0; k < 5; k++) { c[k] = lo[5 * b + k] + hi[5 * b + k]; nn[k] = hi[5 * b + k] - lo[5 * b + k]; }
        civ t0[6][6], tQ[6][6], tW[6][6], tP[6][6];
        for (int k = 0; k < 6; k++) for (int j = 0; j < 6; j++) {
            t0[k][j] = term(k, j, c, nn, pc, ph, 0, v_none);
            tQ[k][j] = term(k, j, c, nn, pc, ph, 0, v_q);
            tW[k][j] = term(k, j, c, nn, pc, ph, 1, v_d);
            tP[k][j] = term(k, j, c, nn, pc, ph, 1, v_all);
        }
        civ S0[6], SQ[6], SW[6], SP[6];
        for (int j = 0; j < 6; j++) {
            S0[j] = SQ[j] = SW[j] = SP[j] = C_ZERO;
            for (int k = 0; k < 6; k++) {
                S0[j] = c_add(S0[j], t0[k][j]); SQ[j] = c_add(SQ[j], tQ[k][j]);
                SW[j] = c_add(SW[j], tW[k][j]); SP[j] = c_add(SP[j], tP[k][j]);
            }
        }
        /* F0, Fq, J, Jq as interval arrays (rows r: j = r+1 for r < 5; constraint rows 5, 6) */
        iv F0[MAXN], Fq[MAXN][MAXP], J[MAXN][MAXN], Jq[MAXP][MAXN][MAXN];
        double Jc[MAXN][MAXN];
        for (int r = 0; r < n; r++) {
            F0[r] = ivp(0.0);
            for (int l = 0; l < MAXP; l++) Fq[r][l] = ivp(0.0);
            for (int s = 0; s < n; s++) { J[r][s] = ivp(0.0); Jc[r][s] = 0.0; for (int l = 0; l < MAXP; l++) Jq[l][r][s] = ivp(0.0); }
        }
        for (int r = 0; r < 5; r++) {
            int j = r + 1;
            F0[r] = iv_sub(c_abs2(S0[j]), ivp(6.0));
            /* Fq: 2 Re(conj(SQ) dSQ/dq_l) */
            for (int l = 0; l < nq; l++) {
                int fl_ = free_[l]; civ d = C_ZERO;
                for (int k = 0; k < 6; k++) d = c_axpy_int(d, -FE[k][j][fl_], c_muli(tQ[k][j]));
                Fq[r][l] = iv_scal(2.0, c_re_conj_mul(SQ[j], d));
            }
            /* J over W (theta, dep vary; q at centre); centre Jacobian Jc for Y */
            for (int kk = 1; kk < 6; kk++) {
                J[r][kk - 1] = iv_scal(2.0, c_re_conj_mul(SW[j], c_muli(tW[kk][j])));
                civ z = c_muli(t0[kk][j]);
                Jc[r][kk - 1] = 2.0 * (midp(S0[j].re) * midp(z.re) + midp(S0[j].im) * midp(z.im));
            }
            for (int i = 0; i < nd; i++) {
                civ d = C_ZERO, dc = C_ZERO;
                for (int k = 0; k < 6; k++) {
                    d = c_axpy_int(d, -FE[k][j][dep[i]], c_muli(tW[k][j]));
                    dc = c_axpy_int(dc, -FE[k][j][dep[i]], c_muli(t0[k][j]));
                }
                J[r][5 + i] = iv_scal(2.0, c_re_conj_mul(SW[j], d));
                Jc[r][5 + i] = 2.0 * (midp(S0[j].re) * midp(dc.re) + midp(S0[j].im) * midp(dc.im));
            }
            /* Jq over the full box */
            for (int l = 0; l < nq; l++) {
                int fl_ = free_[l];
                civ dSq = C_ZERO;
                for (int k = 0; k < 6; k++) dSq = c_axpy_int(dSq, -FE[k][j][fl_], c_muli(tP[k][j]));
                for (int kk = 1; kk < 6; kk++) {
                    /* d/dq_l [2 Re(conj S * i tau_kk)] = 2 Re(conj(dS/dq_l) i tau + conj(S) E_{kk,j,l} tau) */
                    iv a1 = c_re_conj_mul(dSq, c_muli(tP[kk][j]));
                    int e = FE[kk][j][fl_];
                    iv s12 = a1;
                    if (e != 0) s12 = iv_add(a1, c_re_conj_mul(SP[j], c_axpy_int(C_ZERO, e, tP[kk][j])));
                    Jq[l][r][kk - 1] = iv_scal(2.0, s12);
                }
                for (int i = 0; i < nd; i++) {
                    civ dSd = C_ZERO, d2 = C_ZERO;
                    for (int k = 0; k < 6; k++) {
                        dSd = c_axpy_int(dSd, -FE[k][j][dep[i]], c_muli(tP[k][j]));
                        d2 = c_axpy_int(d2, -FE[k][j][fl_] * FE[k][j][dep[i]], tP[k][j]);
                    }
                    iv a1 = c_re_conj_mul(dSq, dSd), a2 = c_re_conj_mul(SP[j], d2);
                    Jq[l][r][5 + i] = iv_scal(2.0, iv_add(a1, a2));
                }
            }
        }
        if (nd) {
            F0[5] = C0.re; F0[6] = C0.im;
            for (int l = 0; l < nq; l++) { Fq[5][l] = dCq[l].re; Fq[6][l] = dCq[l].im; }
            for (int i = 0; i < nd; i++) {
                J[5][5 + i] = dCd[i].re; J[6][5 + i] = dCd[i].im;
                Jc[5][5 + i] = dCd_mid[i][0]; Jc[6][5 + i] = dCd_mid[i][1];
                for (int l = 0; l < nq; l++) { Jq[l][5][5 + i] = d2Cqd[l][i].re; Jq[l][6][5 + i] = d2Cqd[l][i].im; }
            }
        }
        double Y[MAXN][MAXN];
        if (invert(n, Jc, Y) != 0) { excl[b] = 0; if (kd_out) for (int r = 0; r < 2 * MAXN; r++) kd_out[r] = NAN; continue; }
        if (y_out) for (int r = 0; r < n; r++) for (int s2 = 0; s2 < n; s2++) y_out[r * MAXN + s2] = Y[r][s2];
        /* radii of W - w~ */
        double rad[MAXN];
        for (int k = 0; k < 5; k++) rad[k] = angle_up(nn[k]);
        for (int i = 0; i < nd; i++) rad[5 + i] = hq;
        uint8_t ex = 0;
        for (int r = 0; r < n && (!ex || kd_out); r++) {        /* (debug hook: all rows when kd_out is given) */
            /* V_r = (Y F0)_r + sum_l (Y Fq)_{r,l} [-hq, hq] */
            iv V = ivp(0.0);
            for (int s = 0; s < n; s++) V = iv_add(V, iv_scal(Y[r][s], F0[s]));
            for (int l = 0; l < nq; l++) {
                iv yf = ivp(0.0);
                for (int s = 0; s < n; s++) if (!IVZ(Fq[s][l])) yf = iv_add(yf, iv_scal(Y[r][s], Fq[s][l]));
                V = iv_add(V, iv_sym_mul(yf, hq));
            }
            /* Kd_r = -V_r + sum_col E_{r,col} [-rad_col, rad_col],  E = I - Y J - sum_l (Y Jq_l)[-hq,hq] */
            iv Kd = iv_neg(V);
            for (int col = 0; col < n; col++) {
                iv yj = ivp(0.0);
                for (int s = 0; s < n; s++) if (!IVZ(J[s][col])) yj = iv_add(yj, iv_scal(Y[r][s], J[s][col]));
                iv E = iv_sub(ivp(r == col ? 1.0 : 0.0), yj);
                for (int l = 0; l < nq; l++) {
                    iv yjq = ivp(0.0);
                    for (int s = 0; s < n; s++) if (!IVZ(Jq[l][s][col])) yjq = iv_add(yjq, iv_scal(Y[r][s], Jq[l][s][col]));
                    E = iv_add(E, iv_sym_mul(yjq, hq));
                }
                Kd = iv_add(Kd, iv_sym_mul(E, rad[col]));
            }
            if (Kd.lo > rad[r] || Kd.hi < -rad[r]) ex = 1;
            if (kd_out) { kd_out[2 * r] = Kd.lo; kd_out[2 * r + 1] = Kd.hi; }
        }
        excl[b] = ex;
    }
    return 0;
}

/* ---------------------------------------------------------------- mean-value (centred form) exclusion test
 * For each column j: F_j = |S_j|^2 - 6 over the box satisfies (mean value theorem, the box is convex)
 *   F_j(theta, p) in F_j(centre) + sum_k dF_j/dtheta_k(box) [-r_k, r_k] + sum_i dF_j/dp_i(box) [-hq, hq],
 * dF_j/dtheta_k = 2 Re(conj(S_j) i tau_kj),  dF_j/dp_i = 2 Re(conj(S_j) sum_k (-i E_kj,i) tau_kj).
 * Excluded iff this interval misses 0 for some j.  (No constraint is used: valid for every p in the box.) */
int mv_test(int64_t nb, const int64_t *lo, const int64_t *hi, const int64_t *pc, int64_t ph, uint8_t *excl) {
    int v_none[MAXP] = {0, 0, 0, 0}, v_all[MAXP] = {1, 1, 1, 1};
    double hq = angle_up(ph);
    for (int64_t b = 0; b < nb; b++) {
        int64_t c[5], nn[5];
        for (int k = 0; k < 5; k++) { c[k] = lo[5 * b + k] + hi[5 * b + k]; nn[k] = hi[5 * b + k] - lo[5 * b + k]; }
        double rad[5];
        for (int k = 0; k < 5; k++) rad[k] = angle_up(nn[k]);
        uint8_t ex = 0;
        for (int j = 0; j < 6 && !ex; j++) {
            civ t0[6], tP[6], S0 = C_ZERO, SP = C_ZERO;
            for (int k = 0; k < 6; k++) {
                t0[k] = term(k, j, c, nn, pc, ph, 0, v_none); tP[k] = term(k, j, c, nn, pc, ph, 1, v_all);
                S0 = c_add(S0, t0[k]); SP = c_add(SP, tP[k]);
            }
            iv F = iv_sub(c_abs2(S0), ivp(6.0));
            for (int k = 1; k < 6; k++)
                F = iv_add(F, iv_sym_mul(iv_scal(2.0, c_re_conj_mul(SP, c_muli(tP[k]))), rad[k - 1]));
            for (int i = 0; i < FM; i++) {
                civ d = C_ZERO;
                for (int k = 0; k < 6; k++) d = c_axpy_int(d, -FE[k][j][i], c_muli(tP[k]));
                F = iv_add(F, iv_sym_mul(iv_scal(2.0, c_re_conj_mul(SP, d)), hq));
            }
            if (F.lo > 0.0 || F.hi < 0.0) ex = 1;
        }
        excl[b] = ex;
    }
    return 0;
}

int krawczyk_test(int64_t nb, const int64_t *lo, const int64_t *hi, const int64_t *pc, int64_t ph,
                  int nd, const int32_t *dep, uint8_t *excl) {
    return krawczyk_core(nb, lo, hi, pc, ph, nd, dep, excl, 0, 0);
}
/* test hook: for ONE box, the interval vector Kd = K - w~ (all rows) and the preconditioner Y (row-major MAXN) */
int dbg_krawczyk(const int64_t *lo, const int64_t *hi, const int64_t *pc, int64_t ph, int nd, const int32_t *dep,
                 double *kd, double *y) {
    uint8_t ex;
    return krawczyk_core(1, lo, hi, pc, ph, nd, dep, &ex, kd, y);
}

/* ---------------------------------------------------------------- containment of boxes in cyclic hulls (G16) */
int contained(int64_t nb, const int64_t *lo, const int64_t *hi, int64_t nh, const int64_t *hlo, const int64_t *hhi,
              uint8_t *out) {
    const int64_t N16 = 65536;
    for (int64_t b = 0; b < nb; b++) {
        uint8_t in = 0;
        for (int64_t h = 0; h < nh && !in; h++) {
            int ok = 1;
            for (int k = 0; k < 5 && ok; k++) {
                int64_t L = hhi[5 * h + k] - hlo[5 * h + k];
                if (L >= N16) continue;
                int64_t off = ((lo[5 * b + k] - hlo[5 * h + k]) % N16 + N16) % N16;
                if (off + (hi[5 * b + k] - lo[5 * b + k]) > L) ok = 0;
            }
            in = (uint8_t)ok;
        }
        out[b] = in;
    }
    return 0;
}

/* ---------------------------------------------------------------- pair bounds for hulls
 * For v in hull i, w in hull j: <v, w> = sum_k exp(i(phi_k - theta_k)) (k = 1..6, first coordinate 1),
 * phi_k - theta_k in [c_jk - c_ik -+ (n_ik + n_jk)].  |<v,w>|^2 = 0 <=> orthogonal; = 6 <=> unbiased.
 * O[i,j] = 1 iff 0 is possible; B[i,j] = 1 iff 6 is possible.  selfok[i] = 1 iff any two vectors of hull i
 * certainly have |<v,w>|^2 > 6 (then no two distinct vectors of a quartet can share hull i). */
int pair_bounds(int64_t nh, const int64_t *lo, const int64_t *hi, uint8_t *O, uint8_t *B, uint8_t *selfok) {
    for (int64_t i = 0; i < nh; i++) {
        civ S = { { 1.0, 1.0 }, { 0.0, 0.0 } };
        for (int k = 0; k < 5; k++) S = c_add(S, rot(2 * (hi[5 * i + k] - lo[5 * i + k])));
        iv q = c_abs2(S);
        selfok[i] = (uint8_t)(S.re.lo > 0.0 && q.lo > 6.0);
        O[i * nh + i] = 0; B[i * nh + i] = 0;
        for (int64_t j = i + 1; j < nh; j++) {
            civ T = { { 1.0, 1.0 }, { 0.0, 0.0 } };
            for (int k = 0; k < 5; k++) {
                int64_t ci = lo[5 * i + k] + hi[5 * i + k], cj = lo[5 * j + k] + hi[5 * j + k];
                int64_t w = (hi[5 * i + k] - lo[5 * i + k]) + (hi[5 * j + k] - lo[5 * j + k]);
                T = c_add(T, expi(cj - ci, w));
            }
            iv p = c_abs2(T);
            uint8_t o = (uint8_t)(p.lo <= 0.0), bb = (uint8_t)(p.lo <= 6.0 && 6.0 <= p.hi);
            O[i * nh + j] = O[j * nh + i] = o;
            B[i * nh + j] = B[j * nh + i] = bb;
        }
    }
    return 0;
}

/* ---------------------------------------------------------------- Szöllősi parameter selector: (C) and (S)
 * Box of level lev, index idx (4): centre units pc = (2 idx + 1) 2^(16-lev), half-width ph = 2^(16-lev).
 * (C) C(c + d) = C(c) + <grad C(c), d> + R, |R| <= 6 hw^2 (|d_i| <= hw).  Excluded iff for a candidate
 *     direction nv: <nv, -C(c)> - hw sum_i |<nv, D_i>| > 6 hw^2 |nv| (support function of the zonotope).
 * (S) alpha = x + y + z lies within 4 hw (chord <= arc) of alpha(c); excluded iff the disc misses the
 *     closed sector {Im >= 0, Re >= 0, Re - sqrt3 Im >= 0}.  s3lo <= sqrt3 <= s3hi passed in. */
int szollosi_param_test(int64_t nb, const int64_t *idx, int lev, double s3lo, double s3hi, uint8_t *exclC, uint8_t *exclS) {
    int64_t ph = (int64_t)1 << (16 - lev);
    double hw = angle_up(ph);
    double rem = up1(6.0 * up1(hw * hw));
    for (int64_t b = 0; b < nb; b++) {
        int64_t p[4];
        for (int i = 0; i < 4; i++) p[i] = (2 * idx[4 * b + i] + 1) * ph;
        civ x = tab(p[0]), y = tab(p[1]), u = tab(p[2]), v = tab(p[3]);
        civ z = tab(-p[0] - p[1]), w = tab(-p[2] - p[3]);
        civ Cc = c_add(c_add(c_add(x, y), c_add(z, u)), c_add(v, w));
        civ D[4];
        civ mz = c_neg(z), mw = c_neg(w);
        D[0] = c_muli(c_add(x, mz)); D[1] = c_muli(c_add(y, mz)); D[2] = c_muli(c_add(u, mw)); D[3] = c_muli(c_add(v, mw));
        /* candidate directions (numerics only; any direction is valid) */
        double P0 = -midp(Cc.re), P1 = -midp(Cc.im);
        double g[4][2];
        for (int i = 0; i < 4; i++) { g[i][0] = hw * midp(D[i].re); g[i][1] = hw * midp(D[i].im); }
        double cand[25][2]; int nc = 0;
        for (int i = 0; i < 4; i++) { cand[nc][0] = -g[i][1]; cand[nc][1] = g[i][0]; nc++; cand[nc][0] = g[i][1]; cand[nc][1] = -g[i][0]; nc++; }
        for (int s = 0; s < 16; s++) {
            double vx = 0, vy = 0;
            for (int i = 0; i < 4; i++) { double sg = ((s >> i) & 1) ? 1.0 : -1.0; vx += sg * g[i][0]; vy += sg * g[i][1]; }
            cand[nc][0] = P0 - vx; cand[nc][1] = P1 - vy; nc++;
        }
        cand[nc][0] = P0; cand[nc][1] = P1; nc++;
        uint8_t ec = 0;
        for (int a = 0; a < nc && !ec; a++) {
            double n1 = cand[a][0], n2 = cand[a][1];
            if (!(fabs(n1) + fabs(n2) > 0) || !isfinite(n1) || !isfinite(n2)) continue;
            double nn = up1(sqrt(up1(up1(n1 * n1) + up1(n2 * n2))));
            iv dot = iv_add(iv_scal(-n1, Cc.re), iv_scal(-n2, Cc.im));        /* <n, -C(c)> */
            double h = 0.0;
            for (int i = 0; i < 4; i++) h = up1(h + iv_mag(iv_add(iv_scal(n1, D[i].re), iv_scal(n2, D[i].im))));
            double num = dn1(dot.lo - up1(hw * h));
            if (num > up1(nn * rem)) ec = 1;
        }
        exclC[b] = ec;
        /* (S) */
        civ al = c_add(c_add(x, y), z);
        double r = up1(4.0 * hw);
        uint8_t es = 0;
        if (up1(al.im.hi + r) < 0.0) es = 1;                                 /* Im alpha < 0 */
        if (up1(al.re.hi + r) < 0.0) es = 1;                                 /* Re alpha < 0 */
        /* (sqrt3 Im - Re) / 2 > r  on the rectangle al: lower bound of sqrt3 Im - Re */
        double s3im = (al.im.lo >= 0) ? dn1(s3lo * al.im.lo) : dn1(s3hi * al.im.lo);
        double val = dn1(s3im - al.re.hi);
        if (val > up1(2.0 * r)) es = 1;
        exclS[b] = es;
    }
    return 0;
}

/* ---------------------------------------------------------------- test hooks (used only by the test-suite) */
int dbg_arith(int64_t n, const double *a, const double *b, double *out) {
    for (int64_t i = 0; i < n; i++) {
        iv x = ivp(a[i]), y = ivp(b[i]);
        iv s = iv_add(x, y), d = iv_sub(x, y), m = iv_mul(x, y), q = iv_sqr(x);
        double r = sqrt(fabs(a[i]));
        out[10 * i + 0] = s.lo; out[10 * i + 1] = s.hi; out[10 * i + 2] = d.lo; out[10 * i + 3] = d.hi;
        out[10 * i + 4] = m.lo; out[10 * i + 5] = m.hi; out[10 * i + 6] = q.lo; out[10 * i + 7] = q.hi;
        out[10 * i + 8] = dn1(r); out[10 * i + 9] = up1(r);
    }
    return 0;
}
int dbg_expi(int64_t n, const int64_t *c, const int64_t *w, double *out) {
    for (int64_t i = 0; i < n; i++) {
        civ e = expi(c[i], w[i]);
        out[4 * i] = e.re.lo; out[4 * i + 1] = e.re.hi; out[4 * i + 2] = e.im.lo; out[4 * i + 3] = e.im.hi;
    }
    return 0;
}
/* rectangles of S_j (j = 0..5) over the phase box (lo, hi) (G16) and parameter box (pc, ph) (G17) */
int dbg_S(const int64_t *lo, const int64_t *hi, const int64_t *pc, int64_t ph, int vt, int vpmask, double *out) {
    int64_t c[5], n[5];
    int vp[MAXP];
    for (int i = 0; i < MAXP; i++) vp[i] = (vpmask >> i) & 1;
    for (int k = 0; k < 5; k++) { c[k] = lo[k] + hi[k]; n[k] = hi[k] - lo[k]; }
    for (int j = 0; j < 6; j++) {
        civ S = C_ZERO;
        for (int k = 0; k < 6; k++) S = c_add(S, term(k, j, c, n, pc, ph, vt, vp));
        out[4 * j] = S.re.lo; out[4 * j + 1] = S.re.hi; out[4 * j + 2] = S.im.lo; out[4 * j + 3] = S.im.hi;
    }
    return 0;
}
