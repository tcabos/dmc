#include "utils.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include <tuple>
#include <filesystem>
#include <fstream>

void printInfo(const std::string &message)
{
    std::cout << "[INFO] " << message << std::endl;
}

void printError(const std::string &message)
{
    std::cerr << "\033[31m[ERROR] " << message << "\033[0m" << std::endl;
}

void printBlockResults(
    size_t simulTime,
    size_t numWalkers,
    double trialEnergy,
    double blockEnergy,
    double runningAvg,
    double runningError,
    double runningStdDev,
    double blockTimeSeconds)
{

    std::cout << "Tau: " << std::setw(4) << simulTime
              << " | NW: " << std::setw(8) << numWalkers
              << " | Etrial: " << std::fixed << std::setprecision(6) << std::setw(11) << trialEnergy
              << " | Eblock: " << std::fixed << std::setprecision(6) << std::setw(11) << blockEnergy
              << " | Running average: " << std::fixed << std::setprecision(6) << std::setw(11) << runningAvg
              << " | Running error: " << std::fixed << std::setprecision(6) << std::setw(9) << runningError
              << " | Running sigma: " << std::fixed << std::setprecision(6) << std::setw(9) << runningStdDev
              << " | Block time: " << std::fixed << std::setprecision(2) << std::setw(6) << blockTimeSeconds << " s" << std::endl;
}

std::string getCurrentDateTime()
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm localTime = *std::localtime(&now_c);
    std::ostringstream oss;
    oss << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

std::tuple<int, int, int> calculateDuration(std::chrono::high_resolution_clock::time_point start_time,
                                            std::chrono::high_resolution_clock::time_point end_time)
{
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

    int hours = duration / 3600;
    int minutes = (duration % 3600) / 60;
    int seconds = duration % 60;

    return {hours, minutes, seconds};
}

std::string getFormattedTime(int hours, int minutes, int seconds)
{
    std::stringstream ss;
    ss << "Simulation time: "
       << std::setw(2) << std::setfill('0') << hours << ":"
       << std::setw(2) << std::setfill('0') << minutes << ":"
       << std::setw(2) << std::setfill('0') << seconds;
    return ss.str();
}

void writeFinalResultsToFile(DataFile* result_file,
                             double tau, double sim_time, double eq_time, 
                             int block_size, int n_steps, int n_target, 
                             double E_mean_energy, double E_stddev_energy, double E_error_energy,
                             double E_mean_etrial, double E_stddev_etrial, double E_error_etrial,
                             int hours, int minutes, int seconds)
{

    result_file->write_line("tau: " + std::to_string(tau) + " (a.u.)");
    result_file->write_line("simul_time: " + std::to_string(static_cast<size_t>(sim_time)));
    result_file->write_line("equil_time: " + std::to_string(static_cast<size_t>(eq_time)));
    result_file->write_line("block_size: " + std::to_string(block_size));
    result_file->write_line("n_steps: " + std::to_string(n_steps));
    result_file->write_line("n_target: " + std::to_string(n_target));
    result_file->write_line(getFormattedTime(hours, minutes, seconds));
    result_file->write_line("-------------------------------------------------------");
    result_file->write_line("E_mean_energy: " + std::to_string(E_mean_energy));
    result_file->write_line("E_stddev_energy: " + std::to_string(E_stddev_energy));
    result_file->write_line("E_error_energy: " + std::to_string(E_error_energy));
    result_file->write_line("E_final_energy: " + std::to_string(E_mean_energy) + " ± " + std::to_string(E_error_energy));
    result_file->write_line("-------------------------------------------------------");
    result_file->write_line("E_mean_etrial: " + std::to_string(E_mean_etrial));
    result_file->write_line("E_stddev_etrial: " + std::to_string(E_stddev_etrial));
    result_file->write_line("E_error_etrial: " + std::to_string(E_error_etrial));
    result_file->write_line("E_final_etrial: " + std::to_string(E_mean_etrial) + " ± " + std::to_string(E_error_etrial));
}


void snapshot_walker_positions(const std::vector<Walker> &walkers, double time, bool isProduction)
{
    // Určení složky podle přepínače isProduction
    std::string folder = isProduction ? "production" : "equilibration";

    // Vytvoření složky, pokud neexistuje
    if (!std::filesystem::exists(folder))
    {
        std::filesystem::create_directory(folder); // Pokud složka neexistuje, vytvoří ji
    }

    // Vytvoření názvu souboru podle času
    std::ostringstream filename;
    filename << folder << "/walker_snapshot_" << std::setfill('0') << std::setw(2) << static_cast<int>(time) << ".dat";

    // Otevření souboru pro zápis
    std::ofstream file(filename.str(), std::ios::trunc);

    if (!file.is_open())
    {
        throw std::runtime_error("Chyba při otevírání souboru pro zápis.");
    }

    // Zápis pozic walkerů do souboru
    for (size_t i = 0; i < walkers.size(); i++)
    {
        file << walkers[i].e1.x << " " << walkers[i].e1.y << " " << walkers[i].e1.z << " ";
        file << walkers[i].e2.x << " " << walkers[i].e2.y << " " << walkers[i].e2.z << "\n";
    }

    // Soubor je automaticky zavřený při destrukci objektu 'file'
}

