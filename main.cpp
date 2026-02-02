#include "feprog.h"
#include <iostream>

int main(int argc, char *argv[]) {
    std::cout << "TransFE Finite Element Analysis" << std::endl;
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file> [formulation]" << std::endl;
        return 1;
    }

    int formulation = 0;
    if (argc >= 3) {
        formulation = std::atoi(argv[2]);
    }

    try {
        FEProg app;
        app.run_FEA(argv[1], formulation);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
