from ppvm import GeneralizedTableauSimulator
from bloqade import squin

@squin.kernel
def main():
    q = squin.qalloc(2)

    squin.h(q[0])
    squin.cnot(q[0], q[1])

sim = GeneralizedTableauSimulator(2)
task = sim.task(main)
task.run()
print(task.state)
