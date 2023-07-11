#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class CSVParser {
  public:
    explicit CSVParser(std::string _filename, char _delimiter = ',')
        : filename(std::move(_filename)), delimiter(_delimiter) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Failed to open file " << filename << std::endl;
            exit(1);
        }

        std::string line;
        if (std::getline(file, line)) { parse_line(line, header); }

        while (std::getline(file, line)) {
            std::vector<std::string> row;
            parse_line(line, row);
            if (row.size() != header.size()) {
                std::cerr << "Error: Header size (" << header.size() << ") and row size ("
                          << row.size() << ") do not match on line: " << line << std::endl;
                exit(1);
            }
            data.push_back(row);
        }

        file.close();
    }

    std::vector<std::string> get_header() const { return header; }
    size_t nbr_columns() const { return header.size(); }

    bool has_column(const std::string& col) const {
        return find_column_index(col) != std::string::npos;
    }

    std::vector<std::string> get_column(const std::string& col) const {
        std::vector<std::string> column;
        size_t colIndex = find_column_index(col);
        if (colIndex != std::string::npos) {
            for (const auto& row : data) { column.push_back(row[colIndex]); }
        }
        return column;
    }

    size_t find_column_index(const std::string& col) const {
        for (size_t i = 0; i < header.size(); ++i) {
            if (header[i] == col) { return i; }
        }
        return std::string::npos;
    }

    size_t nbr_rows() const { return data.size(); }
    std::vector<std::string> get_row(size_t rowIndex) const { return data.at(rowIndex); }

    std::string get_value(size_t colIndex, size_t rowIndex) const {
        if (colIndex != std::string::npos) { return data.at(rowIndex).at(colIndex); }
        return "";
    }
    double get_value_as_double(size_t colIndex, size_t rowIndex) const {
        string v = get_value(colIndex, rowIndex);
        string v_lower = v;
        std::transform(v_lower.begin(), v_lower.end(), v_lower.begin(), ::tolower);
        if (v.empty() or (v_lower == "nan")) { return std::numeric_limits<double>::quiet_NaN(); }
        return std::stod(v);
    }

  private:
    std::string filename;
    char delimiter;
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> data;

    void parse_line(const std::string& line, std::vector<std::string>& fields) const {
        std::istringstream iss(line);
        std::string field;
        while (std::getline(iss, field, delimiter)) { fields.push_back(field); }
    }
};
