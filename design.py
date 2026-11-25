"""Certified sensor placement on a 31-point grid for a pinned-pinned beam."""

import argparse
import itertools
import json
import math

SCALE = 10**28
GRID = range(1, 32)


def mul(a, b):
    """Outward-rounded product of nonnegative, fixed-point intervals."""
    lo, hi = a[0] * b[0], a[1] * b[1]
    return lo // SCALE, (hi + SCALE - 1) // SCALE


def cosine_table():
    """Bound cos(j*pi/32) using half angles and a recurrence."""
    lo, hi = 0, 0
    for _ in range(4):
        lo = math.isqrt((SCALE + lo) * SCALE // 2)
        hi = math.isqrt((SCALE + hi) * SCALE // 2) + 1
    c = (lo, hi)
    values = [(SCALE, SCALE), c]
    for j in range(2, 16):
        lo, hi = mul(c, values[-1])
        prev_lo, prev_hi = values[-2]
        values.append((2 * lo - prev_hi, 2 * hi - prev_lo))
    values.append((0, 0))
    return values + [(-upper, -lower) for lower, upper in reversed(values[:-1])]


def minors(m):
    cosines = cosine_table()
    differences = {
        (i, j): (cosines[i][0] - cosines[j][1], cosines[i][1] - cosines[j][0])
        for i, j in itertools.combinations(GRID, 2)
    }
    if any(lo <= 0 for lo, _ in differences.values()):
        raise RuntimeError("cosine bounds do not separate all grid points")
    factor = 2 ** (m * (m - 1) // 2)
    result = {}
    for points in itertools.combinations(GRID, m):
        bounds = (factor * SCALE, factor * SCALE)
        for i in points:
            bounds = mul(bounds, cosines[abs(16 - i)])
        for pair in itertools.combinations(points, 2):
            bounds = mul(bounds, differences[pair])
        result[points] = bounds
    return result


class Maximum:
    def __init__(self):
        self.lower = -1
        self.remaining = []
        self.rejected_upper = -1

    def record(self, score, layout, metrics):
        lo, hi = score
        if hi >= self.lower:
            self.remaining.append((hi, layout, metrics))
            self.lower = max(self.lower, lo)
        else:
            self.rejected_upper = max(self.rejected_upper, hi)

    def finalists(self):
        finalists = []
        excluded_upper = self.rejected_upper
        for hi, layout, metrics in self.remaining:
            if hi >= self.lower:
                finalists.append((layout, metrics))
            else:
                excluded_upper = max(excluded_upper, hi)
        return finalists, excluded_upper


def reflection(points):
    return tuple(32 - i for i in reversed(points))


def search(m):
    table = minors(m)
    robust, nominal = Maximum(), Maximum()
    examined = 0
    for p in itertools.combinations(GRID, m + 1):
        deleted = [table[p[:j] + p[j + 1:]] for j in range(m + 1)]
        worst = min(lo for lo, _ in deleted), min(hi for _, hi in deleted)
        gram = sum(lo**2 for lo, _ in deleted), sum(hi**2 for _, hi in deleted)
        metrics = worst, gram
        robust.record(worst, p, metrics)
        nominal.record(gram, p, metrics)
        examined += 1

    if examined != math.comb(len(GRID), m + 1):
        raise RuntimeError("layout count mismatch")

    def summarize(optimum, scale):
        finals, excluded_upper = optimum.finalists()
        positions = [p for p, _ in finals]
        certified = len({min(p, reflection(p)) for p in positions}) == 1
        return {
            "certified_up_to_reflection": certified,
            "score_scale": scale,
            "best_score_lower_scaled": optimum.lower,
            "excluded_layouts_upper_scaled": excluded_upper,
            "layouts": [
            {"positions_over_32": p, "worst_minor": [x / SCALE for x in worst],
             "nominal_gram_det": [x / SCALE**2 for x in gram]}
            for p, (worst, gram) in finals
        ]}

    return {"modes": m, "layouts_examined": examined,
            "robust": summarize(robust, SCALE), "nominal": summarize(nominal, SCALE**2)}


def self_test():
    # Check the product formula against direct determinant expansion.
    for m in range(2, 6):
        points = (3, 9, 16, 24, 29)[:m]
        rows = [[math.sin(k * math.pi * i / 32) for k in range(1, m + 1)] for i in points]
        total = 0
        for perm in itertools.permutations(range(m)):
            inversions = sum(i > j for i, j in itertools.combinations(perm, 2))
            total += (-1)**inversions * math.prod(rows[i][perm[i]] for i in range(m))
        lo, hi = minors(m)[points]
        if not math.isclose(abs(total), (lo + hi) / (2 * SCALE), abs_tol=1e-12):
            raise RuntimeError(f"determinant formula mismatch for m={m}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=int, choices=range(2, 6), nargs="+", default=range(2, 6))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    print(json.dumps([search(m) for m in args.modes], indent=2))


if __name__ == "__main__":
    main()
