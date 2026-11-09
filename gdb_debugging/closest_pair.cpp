#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct Atom { std::string element; double x, y, z; };
struct ClosestPair { std::size_t first, second; double distance; };

std::vector<Atom> read_xyz(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open XYZ file: " + filename);
    int count;
    if (!(input >> count) || count < 2)
        throw std::runtime_error("XYZ must contain at least two atoms");
    std::string line;
    std::getline(input, line);
    std::getline(input, line);
    std::vector<Atom> atoms;
    for (int i = 0; i < count; ++i) {
        Atom atom;
        if (!(input >> atom.element >> atom.x >> atom.y >> atom.z))
            throw std::runtime_error("Invalid XYZ atom record");
        atoms.push_back(atom);
    }
    return atoms;
}

double distance(const Atom& a, const Atom& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

ClosestPair find_closest_pair(const std::vector<Atom>& atoms) {
    const std::size_t n_atoms = atoms.size();
    ClosestPair result{0, 1, std::numeric_limits<double>::infinity()};
    for (std::size_t i = 0; i < n_atoms; ++i) {
        for (std::size_t j = i + 1; j <= n_atoms; ++j) {
            const double separation = distance(atoms.at(i), atoms.at(j));
            if (separation < result.distance) {
                result = {i, j, separation};
            }
        }
    }
    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " coordinates.xyz\n";
        return 1;
    }
    const auto atoms = read_xyz(argv[1]);
    const auto result = find_closest_pair(atoms);
    std::cout << "Closest pair (zero-based): " << result.first << " " << result.second
              << "\nDistance (angstrom): " << std::fixed << std::setprecision(6)
              << result.distance << "\n";
}
