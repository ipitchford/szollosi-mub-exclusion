"""numerical helpers for tests only (not on any certifying path)"""
import numpy as np
def cubic_roots(alpha):
    return np.roots([1.0, -alpha, np.conj(alpha), -1.0])
def discriminant(alpha):
    m = abs(alpha); return m ** 4 + 18 * m ** 2 - 8 * (alpha ** 3).real - 27
def in_region(alpha):
    return discriminant(alpha) <= 0 and discriminant(-alpha) <= 0
def angles_for_alpha(alpha):
    r = cubic_roots(alpha); q = cubic_roots(-alpha)
    r = r[np.argsort(np.angle(r) % (2 * np.pi))]; q = q[np.argsort(np.angle(q) % (2 * np.pi))]
    return np.array([np.angle(t) % (2 * np.pi) for t in (r[0], r[1], q[0], q[1])])
