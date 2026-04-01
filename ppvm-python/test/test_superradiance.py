from ppvm.timeevolve import LadderOp, LindbladOp, solve

from ppvm import PauliSum

n = 4
gamma = 1.0

z = PauliSum.new(n, [f"Z{i}" for i in range(n)])

jump_ops = [LadderOp(i, direction="lower") for i in range(n)]

print(z)

lindblad = LindbladOp(
    jump_ops=jump_ops,
    rates=[gamma] * n,
)


tsteps = 51
tmax = 5.0
tlist = [t / tsteps * tmax for t in range(tsteps)]

print(tlist)

_, values = solve(
    state=z,
    lindblad=lindblad,
    t_span=(0.0, tmax),
    save_at=tlist,
    observable=[f"trace:Z{i}" for i in range(n)],
)

print(values)
