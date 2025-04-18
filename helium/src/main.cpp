#include "statistics.hpp"
#include "data_file.hpp"
#include "utils.hpp"
#include "vec3.hpp"

#include <random>
#include <iostream>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <memory>
#include <sstream>

struct Walker
{
    Vec3 e1;
    Vec3 e2;
    Walker(const Vec3 &e1_ = Vec3(), const Vec3 &e2_ = Vec3()) : e1(e1_), e2(e2_) {}
};

void initialize_walkers(std::vector<Walker> &walkers, size_t target_size)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(0.0, 1.0);

    walkers.resize(target_size);

    for (size_t i = 0; i < walkers.size(); i++)
    {
        walkers[i].e1.x = uniform(gen);
        walkers[i].e1.y = uniform(gen);
        walkers[i].e1.z = uniform(gen);
        walkers[i].e2.x = uniform(gen);
        walkers[i].e2.y = uniform(gen);
        walkers[i].e2.z = uniform(gen);
    }
}

inline double potential(Walker &walker, double rc, double Z)
{
    double r1 = walker.e1.norm();
    double r2 = walker.e2.norm();
    double r12 = walker.e1.distance(walker.e2);

    double V1 = (r1 < rc) ? -(3 * Z) / (2 * rc) : -(Z / r1);
    double V2 = (r2 < rc) ? -(3 * Z) / (2 * rc) : -(Z / r2);

    return V1 + V2 + (1.0 / r12);
}

void simulate(double _tau = 0.01, double _sim_time = 1000, size_t _n_target = 10000, int flag_type_simul = 0)
{
    const size_t n_target = _n_target;
    const double tau = _tau;
    const double sim_time = _sim_time;
    const double eq_time = 30;
    const size_t block_size = static_cast<size_t>(1.0 / tau);
    const size_t n_steps = static_cast<size_t>((sim_time + eq_time) / tau);
    const double Z = 2.0;
    const double D = std::sqrt(tau);
    const size_t m_max = 10;
    const double c_param = 0.1 / tau;
    double e_trial = -2.17;


    std::vector<Walker> walkers;
    std::vector<Walker> temp_walkers;

    std::vector<double> block_local_energies;
    std::vector<double> tmp_local_energies;

    std::vector<double> block_etrial_energies;
    std::vector<double> tmp_etrial_energies;

    std::vector<double> local_energies;

    std::unique_ptr<DataFile> stats_file;
    std::unique_ptr<DataFile> block_energy_file;
    std::unique_ptr<DataFile> block_etrial_file;
    std::unique_ptr<DataFile> result_file;
    std::unique_ptr<DataFile> running_block_file;

    stats_file = std::make_unique<DataFile>("running.dat", 1000);
    block_energy_file = std::make_unique<DataFile>("block_energy.dat", 100);
    block_etrial_file = std::make_unique<DataFile>("block_etrial.dat", 100);
    result_file = std::make_unique<DataFile>("result.dat", 100);
    running_block_file = std::make_unique<DataFile>("running_block.dat", 1);

    walkers.reserve(n_target * 3);
    temp_walkers.reserve(n_target * 3);

    initialize_walkers(walkers, n_target);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (size_t step = 0; step < n_steps; step++)
    {
        double current_time = step * tau - eq_time;
        size_t current_N = walkers.size();

        temp_walkers.clear();
        local_energies.clear();

#pragma omp parallel
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> thread_gauss(0.0, 1.0);
            std::uniform_real_distribution<double> thread_uniform(0.0, 1.0);

            std::vector<Walker> thread_walkers;
            std::vector<double> thread_walkers_local_energies;

#pragma omp for schedule(dynamic, 200) nowait
            for (size_t i = 0; i < current_N; i++)
            {
                Walker old_walker = walkers[i];

                walkers[i].e1.x += D * thread_gauss(gen);
                walkers[i].e1.y += D * thread_gauss(gen);
                walkers[i].e1.z += D * thread_gauss(gen);
                walkers[i].e2.x += D * thread_gauss(gen);
                walkers[i].e2.y += D * thread_gauss(gen);
                walkers[i].e2.z += D * thread_gauss(gen);

                Walker new_walker = walkers[i];

                if (flag_type_simul == 1)
                {
                    double delta_old = old_walker.e1.norm() - old_walker.e2.norm();
                    double delta_new = new_walker.e1.norm() - new_walker.e2.norm();

                    if (delta_old * delta_new < 0.0)
                        continue;
                }

                double A = potential(old_walker, 0, Z);
                double B = potential(new_walker, 0, Z);

                double e_local = (A + B) / 2;

                double w = std::exp(-tau * (e_local - e_trial));
                size_t m = static_cast<size_t>(w + thread_uniform(gen));

                (m > m_max) ? m = m_max : m;

                for (size_t k = 0; k < m; k++)
                {
                    thread_walkers.push_back(new_walker);
                    thread_walkers_local_energies.push_back(B);
                }
            }

#pragma omp critical
            {
                temp_walkers.insert(temp_walkers.end(), thread_walkers.begin(), thread_walkers.end());
                thread_walkers.clear();

                local_energies.insert(local_energies.end(), thread_walkers_local_energies.begin(), thread_walkers_local_energies.end());
                thread_walkers_local_energies.clear();
            }
        }

        std::swap(walkers, temp_walkers);
        current_N = walkers.size();

        if (current_N == 0)
        {
            std::cout << "[ERROR] All walkers have died at time " << step * tau << " au. Exiting.\n";
            break;
        }

        double average_local_energy = Statistics::mean(local_energies);

        e_trial = average_local_energy - c_param * log(((double)current_N / (double)n_target));

        if (step % static_cast<int>(block_size) == 0)
        {
            if (current_time < 0)
            {
                std::ostringstream ss;
                ss << "Time " << current_time << " au (equilibration): N_walkers = " << current_N
                   << ", E = " << average_local_energy << ", E_trial = " << e_trial << "\n";

                std::cout << ss.str();
                running_block_file->write_line(ss.str());
            }
            else
            {
                std::ostringstream ss;
                ss << "Time " << current_time << " au: N_walkers = " << current_N
                   << ", E = " << average_local_energy << ", E_trial = " << e_trial << "\n";

                std::cout << ss.str();
                running_block_file->write_line(ss.str());
            }
        }

        if (current_time >= 0)
        {
            tmp_local_energies.push_back(average_local_energy);
            tmp_etrial_energies.push_back(e_trial);

            if (tmp_local_energies.size() == block_size)
            {
                double block_avg_energy = Statistics::mean(tmp_local_energies);
                double block_avg_Etrial = Statistics::mean(tmp_etrial_energies);

                block_local_energies.push_back(block_avg_energy);
                block_etrial_energies.push_back(block_avg_Etrial);

                block_energy_file->write_numbers(std::vector<double>{
                    block_avg_energy,
                    Statistics::mean(block_local_energies),
                    Statistics::stddev(block_local_energies),
                    Statistics::stderr(block_local_energies)});

                block_etrial_file->write_numbers(std::vector<double>{
                    block_avg_Etrial,
                    Statistics::mean(block_etrial_energies),
                    Statistics::stddev(block_etrial_energies),
                    Statistics::stderr(block_etrial_energies)});

                tmp_local_energies.clear();
                tmp_etrial_energies.clear();
            }
        }

        //stats_file->write_numbers(std::vector<double>{current_time, (double)current_N, average_local_energy, e_trial});
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    int hours, minutes, seconds;
    std::tie(hours, minutes, seconds) = calculateDuration(start_time, end_time);

    if (!block_local_energies.empty())
    {
        double E_mean = Statistics::mean(block_local_energies);
        double Et_mean = Statistics::mean(block_etrial_energies);

        double E_stddev = Statistics::stddev(block_local_energies);
        double Et_stddev = Statistics::stddev(block_etrial_energies);

        double E_error = Statistics::stderr(block_local_energies);
        double Et_error = Statistics::stderr(block_etrial_energies);

        std::cout << "Final Energy Estimate (after equilibration, block averaged):\n";
        std::cout << "<E> = " << E_mean << " ± " << E_error << " , " << E_stddev << std::endl;
        std::cout << "<E_trial> = " << Et_mean << " ± " << Et_error << " , " << Et_stddev << std::endl;
        std::cout << "Simulation time: " << getFormattedTime(hours, minutes, seconds) << std::endl;

        writeFinalResultsToFile(
            result_file.get(), tau, sim_time, eq_time,
            block_size, n_steps, n_target,
            E_mean, E_stddev, E_error,
            Et_mean, Et_stddev, Et_error,
            hours, minutes, seconds);

        result_file->flush();
        stats_file->flush();
        block_energy_file->flush();
        block_etrial_file->flush();
        running_block_file->flush();

        std::cout << "\nSimulation completed successfully. Data written to files.\n";
    }
    else
    {
        std::cout << "\nSimulation terminated early. No statistical data collected.\n";
    }
}

int main(int argc, char const *argv[])
{
    double tau = 0.01;
    double n_target = 10000;
    double sim_time = 100;
    int flag_simul_type = 0;

    if (argc >= 2)
        tau = std::stod(argv[1]);
    if (argc >= 3)
        sim_time = std::stod(argv[2]);
    if (argc >= 4)
        n_target = std::stod(argv[3]);
    if (argc >= 5)
        flag_simul_type = std::stod(argv[4]);

    simulate(tau, sim_time, n_target, flag_simul_type);
    return 0;
}
