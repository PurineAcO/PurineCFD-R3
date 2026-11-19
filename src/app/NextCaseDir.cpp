#include "app/App.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <string>

std::string next_case_dir(const std::string& root, const std::string& tag){
    int next = 1;
    std::error_code ec;
#if chose2D
    const std::string prefix = "2D-";
#else
    const std::string prefix = "3D-";
#endif
    if(std::filesystem::is_directory(root, ec)){
        for(const std::filesystem::directory_entry& e : std::filesystem::directory_iterator(root, ec)){
            const std::string name = e.path().filename().string();
            if(name.compare(0, prefix.size(), prefix) != 0) {continue;}
            const std::size_t dash = name.find('-', prefix.size());
            if(dash == std::string::npos) {continue;}
            const int v = std::atoi(name.substr(prefix.size(), dash - prefix.size()).c_str());
            next = std::max(next, v + 1);
        }
    }
    return root + "/" + prefix + std::to_string(next) + "-" + tag;
}
