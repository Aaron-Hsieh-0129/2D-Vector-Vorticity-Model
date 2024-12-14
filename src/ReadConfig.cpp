#include "ReadConfig.hpp"

std::map<std::string, std::string> vvm_read_config(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, std::string> config;
    std::string line;

    // Read each line from the file
    while (std::getline(file, line)) {
        std::istringstream is_line(line);
        std::string key;
        // Extract the key before '='
        if (std::getline(is_line, key, '=')) {
            std::string value;
            // Extract the value after '='
            if (std::getline(is_line, value)) {
                config[key] = value;  // Store the key-value pair in the map
            }
        }
    }

    return config;  // Return the map with all key-value pairs
}