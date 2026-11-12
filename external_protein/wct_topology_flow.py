import cupy as cp

# ============================================================
# GPU SANITY CHECK
# ============================================================

def gpu_check():
    x = cp.random.rand(5_000_000)
    y = cp.sin(x) * cp.cos(x)
    cp.cuda.runtime.deviceSynchronize()
    print("[OK] CuPy running on GPU")
    print("Device:", cp.cuda.runtime.getDevice())
    print("Mean:", float(cp.mean(y)))


# ============================================================
# GEOMETRY
# ============================================================

def arc_length(points):
    d = cp.diff(points, axis=0, append=points[:1])
    return cp.sum(cp.linalg.norm(d, axis=1))


def curvature(points):
    p_prev = cp.roll(points, 1, axis=0)
    p_next = cp.roll(points, -1, axis=0)

    v1 = points - p_prev
    v2 = p_next - points

    cross = cp.linalg.norm(cp.cross(v1, v2), axis=1)
    denom = cp.linalg.norm(v1, axis=1)**3 + 1e-12
    return cross / denom


def curvature_energy(points):
    k = curvature(points)
    return cp.sum(k * k)


# ============================================================
# SELF-AVOIDANCE (TOPOLOGY ENFORCER)
# ============================================================

def self_avoidance_energy(points, eps=0.3):
    diff = points[:, None, :] - points[None, :, :]
    dist = cp.linalg.norm(diff, axis=-1)
    mask = (dist < eps) & (dist > 1e-6)
    penalty = cp.where(mask, (eps - dist)**2, 0.0)
    return cp.sum(penalty)


# ============================================================
# TOTAL ENERGY
# ============================================================

def total_energy(points, sa_weight=50.0):
    return curvature_energy(points) + sa_weight * self_avoidance_energy(points)


# ============================================================
# FINITE-DIFFERENCE GRADIENT (GPU SAFE)
# ============================================================

def energy_gradient(points, sa_weight=50.0, h=1e-3):
    grad = cp.zeros_like(points)

    base_E = total_energy(points, sa_weight)

    for axis in range(3):
        shift = cp.zeros_like(points)
        shift[:, axis] = h

        e_plus = total_energy(points + shift, sa_weight)
        e_minus = total_energy(points - shift, sa_weight)

        grad[:, axis] = (e_plus - e_minus) / (2 * h)

    return grad


# ============================================================
# WCT CURVATURE FLOW
# ============================================================

def evolve(points, steps=500, dt=1e-4, sa_weight=50.0):
    for step in range(steps):
        grad = energy_gradient(points, sa_weight)
        points -= dt * grad

        if step % 25 == 0:
            E = total_energy(points, sa_weight)
            k_mean = cp.mean(curvature(points))
            L = arc_length(points)
            print(
                f"step {step:4d} | "
                f"E={float(E):.3f} | "
                f"<κ>={float(k_mean):.4f} | "
                f"L={float(L):.3f}"
            )

    return points


# ============================================================
# TOPOLOGY SEEDS
# ============================================================

def circle(N=400, R=4.0):
    t = cp.linspace(0, 2 * cp.pi, N, endpoint=False)
    x = R * cp.cos(t)
    y = R * cp.sin(t)
    z = cp.zeros_like(t)
    return cp.stack([x, y, z], axis=1)


def trefoil(N=400, scale=4.0):
    t = cp.linspace(0, 2 * cp.pi, N, endpoint=False)
    x = cp.sin(t) + 2 * cp.sin(2 * t)
    y = cp.cos(t) - 2 * cp.cos(2 * t)
    z = -cp.sin(3 * t)
    return scale * cp.stack([x, y, z], axis=1)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("\n=== WCT GPU TOPOLOGY FLOW ===\n")

    gpu_check()

    print("\n--- UNKNOT TEST (should collapse) ---")
    unknot = circle()
    print("Initial energy:", float(total_energy(unknot)))
    unknot = evolve(unknot, steps=200)
    print("Final energy:", float(total_energy(unknot)))

    print("\n--- TREFOIL TEST (should NOT unknot) ---")
    knot = trefoil()
    print("Initial energy:", float(total_energy(knot)))
    knot = evolve(knot, steps=300)
    print("Final energy:", float(total_energy(knot)))

    print("\n[Done]")
