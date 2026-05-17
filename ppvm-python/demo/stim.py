import time
from pathlib import Path

from ppvm import StimProgram, sample_stim

file_path = str(Path(__file__).parent / "msd.stim")

n_shots = 100

start = time.time()
prog = StimProgram.from_file(file_path)
shot_results = sample_stim(prog, n_qubits=85, num_shots=n_shots, seed=0)

runtime = time.time() - start

print(f"Overall runtime for {n_shots} shots of the 85 qubit MSD circuit: {runtime} s")

# can also run from a string
from ppvm import GeneralizedTableau

tab = GeneralizedTableau(2)
prog = StimProgram.parse("""
H 0
CX 0 1
M 0 1
""")
results = tab.run(prog)
print(f"Bell state results: {results}")
