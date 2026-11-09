// Profile minimum interatomic distances for reproducible synthetic configurations.
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct Point { double x, y, z; };
using Structure = std::vector<Point>;
using Matrix = std::vector<std::vector<double>>;

std::vector<Structure> make_structures(int count, int atoms) {
    std::mt19937 generator(274);
    std::uniform_real_distribution<double> jitter(-0.1, 0.1);
    const int side = static_cast<int>(std::ceil(std::cbrt(atoms)));
    std::vector<Structure> structures;
    for (int k = 0; k < count; ++k) {
        Structure coordinates;
        for (int i = 0; i < atoms; ++i) {
            coordinates.push_back({2.0*(i % side) + jitter(generator),
                2.0*((i / side) % side) + jitter(generator),
                2.0*(i / (side*side)) + jitter(generator)});
        }
        structures.push_back(coordinates);
    }
    return structures;
}

Matrix pairwise_distances(const Structure& coordinates) {
    const std::size_t n = coordinates.size();
    Matrix matrix(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double dx = coordinates[i].x - coordinates[j].x;
            const double dy = coordinates[i].y - coordinates[j].y;
            const double dz = coordinates[i].z - coordinates[j].z;
            matrix[i][j] = matrix[j][i] = std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    return matrix;
}

double minimum_distance(const Structure& coordinates) {
    double closest = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        const Matrix distances = pairwise_distances(coordinates);
        for (std::size_t j = i + 1; j < coordinates.size(); ++j) {
            if (distances[i][j] < closest) closest = distances[i][j];
        }
    }
    return closest;
}

std::vector<double> analyze(const std::vector<Structure>& structures) {
    std::vector<double> results;
    for (const Structure& coordinates : structures)
        results.push_back(minimum_distance(coordinates));
    return results;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " atoms structures repeats output.txt\n";
        return 1;
    }
    try {
        const int atoms = std::stoi(argv[1]);
        const int count = std::stoi(argv[2]);
        const int repeats = std::stoi(argv[3]);
        if (atoms < 2 || count < 1 || repeats < 1)
            throw std::runtime_error("Need atoms >= 2, structures >= 1, repeats >= 1");
        const auto structures = make_structures(count, atoms);
        std::vector<double> results;
        double checksum = 0.0;
        const auto start = std::chrono::steady_clock::now();
        for (int r = 0; r < repeats; ++r) {
            results = analyze(structures);
            for (double value : results) checksum += value;
        }
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        std::ofstream output(argv[4]);
        if (!output) throw std::runtime_error("Cannot open output file");
        output << std::setprecision(17);
        for (double value : results) output << value << '\n';
        std::cout << std::setprecision(10) << "Analysis seconds: " << elapsed
                  << "\nChecksum: " << checksum << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
