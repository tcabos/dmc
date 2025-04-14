#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <vector>

class Statistics {
public:
    static double mean(const std::vector<double>& data);
    static double variance(const std::vector<double>& data);
    static double stddev(const std::vector<double>& data);
    static double stderr(const std::vector<double>& data); 
};

#endif // STATISTICS_HPP