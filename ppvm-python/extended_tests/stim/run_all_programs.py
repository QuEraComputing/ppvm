import os
import time
from pathlib import Path

from ppvm import GeneralizedTableau, StimProgram

folder = str(Path(__file__).parent / "stim-programs/")

n_qubits = {}

for file in os.listdir(folder):
    if not file.endswith(".stim"):
        continue

    max_qubit_address = -1
    file_path = os.path.join(folder, file)
    with open(file_path) as contents:
        for line in contents.readlines():
            parts = line.split(" ")
            for part in parts:
                if not part.isdigit():
                    continue

                num = int(part)
                if num > max_qubit_address:
                    max_qubit_address = num

    if max_qubit_address < 0:
        raise ValueError(f"Couldn't find qubit number for {file}")

    n_qubits[file] = max_qubit_address + 1


for file in os.listdir(folder):
    if not file.endswith(".stim"):
        continue

    file_path = os.path.join(folder, file)

    print("=" * 30)
    print(f"Running {file}")

    tab = GeneralizedTableau(n_qubits[file])

    start = time.time()
    tab.run(StimProgram.from_file(file_path))

    print(f"Finished running {file}")
    print(f"Runtime (n_qubits = {n_qubits[file]}): {(time.time() - start) * 1e3} ms")

    print("=" * 30)
