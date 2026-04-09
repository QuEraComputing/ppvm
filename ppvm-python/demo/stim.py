import time
from pathlib import Path

from ppvm import GeneralizedTableau

file_path = str(Path(__file__).parent / "msd.stim")

tab = GeneralizedTableau(85)

shots = 100

start = time.time()
shot_results = []

for _ in range(shots):
    # Create a copy since the tableau is mutated when running the circuit
    tab_shot = tab.fork()
    results = tab_shot.run_stim_file(file_path)
    shot_results.append(results)

runtime = time.time() - start

print(f"Overall runtime for {shots} shots of the 85 qubit MSD circuit: {runtime} s")

# can also run strings
tab = GeneralizedTableau(2)

stim_str = """
H 0
CX 0 1
M 0 1
"""

results = tab.run_stim_string(stim_str)
print(f"Bell state results: {results}")
