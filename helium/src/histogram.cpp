#include "histogram.hpp"
#include <fstream>
#include <cmath>

Histogram::Histogram(double bin_size, double max_range) : bin_size_(bin_size) {
    num_bins_ = static_cast<size_t>(std::ceil(max_range / bin_size));
    bins_.resize(num_bins_, 0);
}

void Histogram::add(double r) {
    size_t bin = static_cast<size_t>(std::floor(r / bin_size_));
    if (bin < num_bins_) ++bins_[bin];
}

void Histogram::writeToFile(const std::string& filename) const {
    std::ofstream fout(filename);
    if (!fout) throw std::runtime_error("Unable to open histogram file");
    for (size_t i = 1; i <= num_bins_; ++i) {
        double bin_center = (i + 0.5) * bin_size_;
        //double bin_center = i * bin_size_;
        fout << bin_center << " " << bins_[i] << "\n";
    }
    fout.close();
}