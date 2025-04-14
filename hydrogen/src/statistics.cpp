#include "statistics.hpp"
#include <numeric>
#include <cmath>
#include <stdexcept>

double Statistics::mean(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Data is empty.");

    double mean = 0.0;
    size_t count = 0;

    for (double x : data) {
        ++count;
        mean += (x - mean) / count;
    }

    return mean;
}

double Statistics::variance(const std::vector<double>& data) {
    if (data.size() < 2) return 0;

    double mean = 0.0;
    double M2 = 0.0;
    size_t count = 0;

    for (double x : data) {
        ++count;
        double delta = x - mean;
        mean += delta / count;
        M2 += delta * (x - mean);
    }

    return M2 / (count - 1); // sample variance
}

double Statistics::stddev(const std::vector<double>& data) {
    return std::sqrt(variance(data));
}

double Statistics::stderr(const std::vector<double>& data) {
    return stddev(data) / std::sqrt(static_cast<double>(data.size()));
}