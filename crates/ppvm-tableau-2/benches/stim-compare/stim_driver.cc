// Times Stim's TableauSimulator on a circuit the way ppvm-tableau-2's run_stim
// runs it: noise removed, a fresh simulator per shot. Builds against any Stim
// checkout's libstim; see README.md. Usage: stim_driver <circuit.stim> [shots].
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

#include "stim.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: stim_driver <circuit.stim> [shots]\n";
        return 1;
    }
    std::ifstream in(argv[1]);
    std::stringstream text;
    text << in.rdbuf();
    stim::Circuit circuit = stim::Circuit(text.str()).without_noise();
    size_t shots = argc > 2 ? std::stoul(argv[2]) : 1000;
    size_t n = circuit.count_qubits();

    std::mt19937_64 rng(0x9e3779b97f4a7c15ULL);
    bool parity = false;
    auto start = std::chrono::steady_clock::now();
    for (size_t s = 0; s < shots; s++) {
        stim::TableauSimulator<stim::MAX_BITWORD_WIDTH> sim(std::move(rng), n);
        sim.safe_do_circuit(circuit);
        if (!sim.measurement_record.storage.empty()) {
            parity ^= sim.measurement_record.storage.back();
        }
        rng = std::move(sim.rng);
    }
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << argv[1] << ": " << n << " qubits, " << shots << " shots, word width "
              << stim::MAX_BITWORD_WIDTH << "\n"
              << elapsed / shots * 1e6 << " us/shot (parity " << parity << ")\n";
}
