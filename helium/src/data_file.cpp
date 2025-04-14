#include "data_file.hpp"
#include <fstream>
#include <iomanip>

DataFile::DataFile(const std::string &name, size_t bufferlimit) : filename(name), bufferlimit(bufferlimit), lineCount(0)
{
    file.open(filename, std::ios::trunc);
    if (!file)
        throw std::runtime_error("Unable to open file " + filename);
}

DataFile::~DataFile()
{
    flush();
    file.close();
}

void DataFile::write_line(const std::string &line)
{
    buffer << line << "\n";
    lineCount++;
    if (lineCount >= bufferlimit)
        flush();
}

void DataFile::write_numbers(const std::vector<double> &data)
{
    for (size_t i = 0; i < data.size(); i++)
        buffer << data[i] << " ";
    buffer << "\n";
    lineCount++;
    if (lineCount >= bufferlimit)
        flush();
}

void DataFile::write_numbers(const std::vector<size_t> &data)
{
    for (size_t i = 0; i < data.size(); i++)
        buffer << data[i] << " ";
    buffer << "\n";
    lineCount++;
    if (lineCount >= bufferlimit)
        flush();
}

void DataFile::write_number(const double &number)
{
    buffer << std::fixed << std::setprecision(6) << number << "\n";
    lineCount++;
    if (lineCount >= bufferlimit)
    flush();
}

void DataFile::write_number(const size_t &number)
{
    buffer << number << "\n";
    lineCount++;
    if (lineCount >= bufferlimit)
    flush();
}

void write_number(const size_t &number);


void DataFile::flush()
{
    std::string buffer_content = buffer.str();
    if (!buffer_content.empty())
    {
        file << buffer_content;
        file.flush();
    }
    buffer.str("");
    lineCount = 0;
}
