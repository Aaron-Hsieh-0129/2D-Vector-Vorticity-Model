#include "ReadConfig.hpp"

std::map<std::string, std::string> vvm_read_config(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, std::string> config;
    std::string line;

    // Read each line from the file
    while (std::getline(file, line)) {
        // Ignore everything after a '#' character
        auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // Remove leading and trailing whitespace (optional)
        line.erase(0, line.find_first_not_of(" \t")); // Leading
        line.erase(line.find_last_not_of(" \t") + 1); // Trailing

        if (line.empty()) {
            continue; // Skip empty lines
        }

        std::istringstream is_line(line);
        std::string key;

        // Extract the key before '='
        if (std::getline(is_line, key, '=')) {
            std::string value;

            // Remove leading/trailing whitespace from the key
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);

            // Extract the value after '='
            if (std::getline(is_line, value)) {
                // Remove leading/trailing whitespace from the value
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);

                // Store the key-value pair in the map
                config[key] = value;
            }
        }
    }

    return config; // Return the map with all key-value pairs
}