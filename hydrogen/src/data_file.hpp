#ifndef DATA_FILE_HPP
#define DATA_FILE_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

class DataFile {
private:
    std::ofstream file;
    std::string filename;
    std::ostringstream buffer;
    size_t bufferlimit;
    size_t lineCount;

public:
    DataFile(const std::string& name, size_t bufferlimit);
    ~DataFile();
    void write_line(const std::string& line);
    void write_numbers(const std::vector<double> &data);
    void write_numbers(const std::vector<size_t> &data);
    void write_number(const double &number);
    void write_number(const size_t &number);
    void flush();
};

#endif // DATA_FILE_HPP