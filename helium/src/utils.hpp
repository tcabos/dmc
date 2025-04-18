#ifndef UTILS_HPP
#define UTILS_HPP

#include "data_file.hpp"
#include "walker.hpp"
#include <string>
#include <tuple>
#include <chrono>

void printInfo(const std::string &message);
void printError(const std::string &message);
void printBlockResults(
    size_t simulTime,
    size_t numWalkers,
    double trialEnergy,
    double blockEnergy,
    double runningAvg,
    double runningError,
    double runningStdDev,
    double blockTimeSeconds);
std::string getCurrentDateTime();
std::string getFormattedTime(int hours, int minutes, int seconds);
void writeFinalResultsToFile(DataFile *result_file,
                             double tau, double sim_time, double eq_time,
                             int block_size, int n_steps, int n_target,
                             double E_mean_energy, double E_stddev_energy, double E_error_energy,
                             double E_mean_etrial, double E_stddev_etrial, double E_error_etrial,
                             int hours, int minutes, int seconds);
std::tuple<int, int, int> calculateDuration(
    std::chrono::high_resolution_clock::time_point start,
    std::chrono::high_resolution_clock::time_point end);
void snapshot_walker_positions(const std::vector<Walker> &walkers, double time, bool isProduction);
#endif // UTILS_HPP