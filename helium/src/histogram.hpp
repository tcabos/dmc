#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include <vector>
#include <string>

class Histogram {
private:
    std::vector<size_t> bins_;
    size_t num_bins_;
    double bin_size_;

public:
    Histogram(double bin_size = 1.0, double max_range = 10.0);
    void add(double r);
    void writeToFile(const std::string& filename) const;
};

#endif // HISTOGRAM_HPP