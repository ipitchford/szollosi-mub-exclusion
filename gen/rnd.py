"""
rnd.py -- floating-point model and outward-rounding helpers of the GENERATOR (untrusted witness producer;
the proof is carried by the independent checker in ../chk, which uses a different arithmetic model).

Floating-point model: IEEE-754 binary64, round-to-nearest, unit roundoff u = 2^-53.  U = 2u = 2^-52.

Two devices make every computed radius a rigorous upper bound (referee comment 4.1):

(1) up(x): per-primitive outward correction.
    LEMMA R.  Let E be an expression built from nonnegative real inputs a_1..a_n by at most K <= 256 operations
    + and x (including matrix products of nonnegative matrices, counted as their inner dimension + 1), with
    multiplication depth <= 4.  Let a_i-hat be floats with a_i <= a_i-hat (1-u)^-8 (i.e. each input is either
    exact or an upper bound, or a floating value within 8 roundings of one), all <= 2^64.  If E-hat is the
    value of E computed in binary64 RN (any evaluation order, BLAS allowed: for nonnegative data a computed
    dot product satisfies fl >= (1-gamma_p) * exact, Higham 2002, Thm 3.5), then
          E(a) <= up(E-hat),   up(x) := fl(fl(x * (1 + 2^-40)) + 2^-800)   (and +inf if x > 2^64 or x is NaN).
    Proof sketch (full proof in the paper, Appendix A): every + or x of nonnegatives loses at most a factor
    (1-u) (additions are exact in the subnormal range; each multiplication may in addition lose an absolute
    2^-1075 to underflow, whose influence on the result is at most 2^-1075 (2^64)^4 = 2^-819 per operation);
    hence E(a) <= E-hat (1-u)^-(K+8) + 256 * 2^-819 <= E-hat (1 + 2^-41) + 2^-811, and
    fl(fl(x(1+2^-40)) + 2^-800) >= (x (1+2^-40)(1-u) + 2^-800)(1-u) exceeds this.
    Every radius returned by a primitive in this package is passed through up().

(2) nu(x) / nd(x): one-ulp outward steps (numpy.nextafter) after a SINGLE correctly rounded operation whose
    operands are exact floats; used for the few non-monotone quantities (differences, lower bounds).
    In RN the exact result lies within half a spacing of the computed one on either side, and the neighbour
    in each direction is at least that far away, so nd(fl(a op b)) <= a op b <= nu(fl(a op b)).

All exclusion decisions are strict comparisons of the form  lower_bound > upper_bound  with both sides
produced by (1)-(2); NaN/inf make every such comparison False, i.e. never lead to an exclusion.
"""
import numpy as np

U = 2.0 ** -52                      # safe unit roundoff (2u)
TAB_ERR = 2.0 ** -52                # |T[m] - exp(2 pi i m / 2^17)| <= sqrt(2) 2^-53 < 2^-52 (Arb-verified)
UPF = 1.0 + 2.0 ** -40
TINY = 2.0 ** -800
BIG = 2.0 ** 64
PI_UP = float(np.nextafter(np.pi, np.inf))       # np.pi < pi < PI_UP  (np.pi = RN(pi) is below pi)
PI_DN = float(np.pi)
assert PI_DN < PI_UP
SQ6 = float(np.sqrt(6.0))                          # RN(sqrt 6); |SQ6 - sqrt6| <= u sqrt6
SQ6_UP = float(np.nextafter(SQ6, np.inf)); SQ6_DN = float(np.nextafter(SQ6, -np.inf))
INV_SQ6_UP = float(np.nextafter(1.0 / SQ6_DN, np.inf))   # >= 1/sqrt6
INV_SQ6_DN = float(np.nextafter(1.0 / SQ6_UP, -np.inf))  # <= 1/sqrt6

def up(x):
    """rigorous upper bound of a monotone nonnegative expression computed in RN (Lemma R)"""
    x = np.asarray(x, dtype=float)
    with np.errstate(over='ignore', invalid='ignore'):
        y = x * UPF + TINY
    return np.where(x <= BIG, y, np.inf)

def nu(x):
    return np.nextafter(x, np.inf)

def nd(x):
    return np.nextafter(x, -np.inf)

def angle_up(n_grid, grid_log2):
    """rigorous upper bound of n * 2 pi / 2^grid_log2 (n >= 0 integer array): the exact value is n*pi/2^(g-1);
    PI_UP / 2^(g-1) is exact (power-of-two scaling), and one product is corrected by nu()."""
    return nu(np.asarray(n_grid, dtype=float) * (PI_UP / 2.0 ** (grid_log2 - 1)))
