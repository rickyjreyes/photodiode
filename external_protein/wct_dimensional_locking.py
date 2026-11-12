import cupy as cp
import numpy as np

# ============================================================
# CONFIG
# ============================================================

DIMENSIONS = [3, 4]
TRIALS = 512
N = 64
STEPS = 1000
DT = 5e-4

TOL = 0.25
TWO_PI = 2 * cp.pi

cp.random.seed(0)

# ============================================================
# RANDOM CLOSED LOOP
# ============================================================

def random_loop(N, dim):
    x = cp.cumsum(cp.random.randn(N, dim), axis=0)
    x -= cp.mean(x, axis=0)
    return x

# ============================================================
# CURVATURE FLOW (DIM-INDEPENDENT)
# ============================================================

def relax(x, steps=STEPS):
    for _ in range(steps):
        dx = cp.roll(x, -1, 0) - x
        ddx = cp.roll(x, -1, 0) - 2*x + cp.roll(x, 1, 0)

        grad = dx * ddx
        x -= DT * grad
        x -= cp.mean(x, axis=0)

    return x

# ============================================================
# FRAME TRANSPORT (GENERALIZED)
# ============================================================

def transport_angle(x):
    """
    Measures accumulated rotation of tangent frames.
    Works in any dimension.
    """
    dx = cp.roll(x, -1, 0) - x
    t = dx / (cp.linalg.norm(dx, axis=1, keepdims=True) + 1e-8)

    dots = cp.sum(t * cp.roll(t, -1, 0), axis=1)
    dots = cp.clip(dots, -1.0, 1.0)

    angles = cp.arccos(dots)
    return cp.sum(angles)

def transport_closes(x):
    phi = transport_angle(x)
    return cp.abs(cp.mod(phi, TWO_PI)) < TOL

# ============================================================
# EXPERIMENT
# ============================================================

def run_experiment(dim):
    survivors = 0

    for _ in range(TRIALS):
        x = random_loop(N, dim)
        x = relax(x)

        if transport_closes(x):
            survivors += 1

    return survivors

# ============================================================
# MAIN
# ============================================================

def main():
    print("\n=== WCT DIMENSIONAL LOCKING TEST ===\n")

    for dim in DIMENSIONS:
        surv = run_experiment(dim)
        print(f"Dimension {dim}D:")
        print(f"  Survivors        : {surv}")
        print(f"  Survival fraction: {surv / TRIALS:.4f}\n")

    print("Interpretation:")
    print("- Finite survival in 3D → transport closure allowed")
    print("- Vanishing survival in ≥4D → dimensional locking")
    print("- If 4D survives comparably → WCT falsified\n")

if __name__ == "__main__":
    main()
