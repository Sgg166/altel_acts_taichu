#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <regex>
#include <numeric>
#include <tuple>
#include "myrapidjson.h"
 
class ParameterOptimizer {
private:
    std::string outputDir;
    std::string dataFile;
    std::string geometryFile;
    std::string configFile;
    std::string rootFile;
    int targetId;
    int eventSkip;
    int eventMax;
    double convergenceThreshold;
    int totalIterations = 0;
    
    // Store the result of each iteration
    struct IterationResult {
        int iteration_num;
        double u_res1, u_res2, u_res3;      
        double v_res1, v_res2, v_res3;      
        double u_err_mean, v_err_mean;      
        double u_pixel1, u_pixel2, u_pixel3;
        double v_pixel1, v_pixel2, v_pixel3; 
    };
    
    std::vector<IterationResult> iterationResults; // Store all iteration results
    
    struct ResolutionParams {
        double u_size1, u_size2, u_sizeOther;
        double v_size1, v_size2, v_sizeOther;
 
        ResolutionParams() : u_size1(6), u_size2(6), u_sizeOther(6),
                          v_size1(6), v_size2(6), v_sizeOther(6) {}
 
        void loadFromJson(const std::string& filename) {
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Cannot open config file: " << filename << std::endl;
                return;
            }
 
            std::string content((std::istreambuf_iterator<char>(file)),std::istreambuf_iterator<char>());
            file.close();
 
            JsonDocument doc;
            if (doc.Parse(content.c_str()).HasParseError()) {
                std::cerr << "JSON parse error in file: " << filename << std::endl;
                return;
            }
 
            std::cout << "Loading config from: " << filename << std::endl;
 
            if (doc.HasMember("resolutions")) {
                const auto& res = doc["resolutions"];
                if (res.HasMember("u_direction")) {
                    const auto& u = res["u_direction"];
                    if (u.HasMember("size_1")) {
                        u_size1 = std::stod(u["size_1"].GetString());
                        std::cout << "u_size1 = " << u_size1 << std::endl;
                    }
                    if (u.HasMember("size_2")) {
                        u_size2 = std::stod(u["size_2"].GetString());
                        std::cout << "u_size2 = " << u_size2 << std::endl;
                    }
                    if (u.HasMember("size_other")) {
                        u_sizeOther = std::stod(u["size_other"].GetString());
                        std::cout << "u_sizeOther = " << u_sizeOther << std::endl;
                    }
                }
                if (res.HasMember("v_direction")) {
                    const auto& v = res["v_direction"];
                    if (v.HasMember("size_1")) {
                        v_size1 = std::stod(v["size_1"].GetString());
                        std::cout << "v_size1 = " << v_size1 << std::endl;
                    }
                    if (v.HasMember("size_2")) {
                        v_size2 = std::stod(v["size_2"].GetString());
                        std::cout << "v_size2 = " << v_size2 << std::endl;
                    }
                    if (v.HasMember("size_other")) {
                        v_sizeOther = std::stod(v["size_other"].GetString());
                        std::cout << "v_sizeOther = " << v_sizeOther << std::endl;
                    }
                }
            }
        }
 
        void saveToJson(const std::string& filename) const {
            JsonDocument doc;
            doc.SetObject();
            JsonAllocator& allocator = doc.GetAllocator();
 
            JsonValue resolutions(rapidjson::kObjectType);
            JsonValue u_direction(rapidjson::kObjectType);
            JsonValue v_direction(rapidjson::kObjectType);
 
            u_direction.AddMember("size_1", JsonValue().SetString(std::to_string(u_size1).c_str(), allocator), allocator);
            u_direction.AddMember("size_2", JsonValue().SetString(std::to_string(u_size2).c_str(), allocator), allocator);
            u_direction.AddMember("size_other", JsonValue().SetString(std::to_string(u_sizeOther).c_str(), allocator), allocator);
 
            v_direction.AddMember("size_1", JsonValue().SetString(std::to_string(v_size1).c_str(), allocator), allocator);
            v_direction.AddMember("size_2", JsonValue().SetString(std::to_string(v_size2).c_str(), allocator), allocator);
            v_direction.AddMember("size_other", JsonValue().SetString(std::to_string(v_sizeOther).c_str(), allocator), allocator);
 
            resolutions.AddMember("u_direction", u_direction, allocator);
            resolutions.AddMember("v_direction", v_direction, allocator);
            doc.AddMember("resolutions", resolutions, allocator);
 
            std::ofstream file(filename);
            rapidjson::StringBuffer buffer;
            rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
            doc.Accept(writer);
            file << buffer.GetString();
        }
 
        bool isConverged(const ResolutionParams& other, double threshold) const {
            double u_diff1 = std::abs(u_size1 - other.u_size1);
            double u_diff2 = std::abs(u_size2 - other.u_size2);
            double u_diff3 = std::abs(u_sizeOther - other.u_sizeOther);
            double v_diff1 = std::abs(v_size1 - other.v_size1);
            double v_diff2 = std::abs(v_size2 - other.v_size2);
            double v_diff3 = std::abs(v_sizeOther - other.v_sizeOther);
 
            return (u_diff1 < threshold && u_diff2 < threshold && u_diff3 < threshold &&
                    v_diff1 < threshold && v_diff2 < threshold && v_diff3 < threshold);
        }
    };
 
public:
    ParameterOptimizer(const std::vector<std::string>& dataPaths, const std::string& geomPath,
                      const std::string& configPath, const std::string& outputPath,
                      int target, int maxEvents = 100000, double convThreshold = 0.000001)
        : geometryFile(geomPath), configFile(configPath),
          rootFile(outputPath), targetId(target), eventSkip(0), eventMax(maxEvents),
          convergenceThreshold(convThreshold) {

        //Merge all data files into a single string, separated by Spaces
        for (size_t i = 0; i < dataPaths.size(); i++) {
            if (i > 0) dataFile += " ";
            dataFile += dataPaths[i];
        }
        outputDir = std::filesystem::path(configPath).parent_path().string() + "/optimization_results";
        std::filesystem::create_directories(outputDir);
    }
    int getTotalIterations() const {
        return totalIterations;
    }
 
   //Obtain the σres values of the three cluster sizes
    std::tuple<double, double, double, double, double, double> getResidualValues() {
        std::ostringstream cmd;

        //这是分单，双，多像素计算残差
        cmd << "root -l -q -b 'Match.C(\"" << rootFile << "\")'";
 
        std::cout << "Running Match.C with file: " << rootFile << std::endl;
 
        FILE* pipe = popen(cmd.str().c_str(), "r");
        std::string output;
        if (pipe) {
            char buffer[256];
            while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                output += buffer;
            }
            pclose(pipe);
        }
 
        std::cout << "Match.C output:" << std::endl << output << std::endl;
 
        // Match the σres values of the three cluster sizes
        std::regex u1_pattern(R"(SigmaU1\s*=\s*([\d.]+)\s*um)");  // Single pixel
        std::regex u2_pattern(R"(SigmaU2\s*=\s*([\d.]+)\s*um)");  // Two pixels
        std::regex u3_pattern(R"(SigmaU3\s*=\s*([\d.]+)\s*um)");  // Multi-pixel
 
        std::regex v1_pattern(R"(SigmaV1\s*=\s*([\d.]+)\s*um)");  
        std::regex v2_pattern(R"(SigmaV2\s*=\s*([\d.]+)\s*um)");  
        std::regex v3_pattern(R"(SigmaV3\s*=\s*([\d.]+)\s*um)");  
 
        std::smatch match;
        double u_res1 = 0.0, u_res2 = 0.0, u_res3 = 0.0;
        double v_res1 = 0.0, v_res2 = 0.0, v_res3 = 0.0;
 
        if (std::regex_search(output, match, u1_pattern)) {
            u_res1 = std::stod(match[1].str());
            std::cout << "Found σres_U1: " << u_res1 << " um" << std::endl;
        }
 
        if (std::regex_search(output, match, u2_pattern)) {
            u_res2 = std::stod(match[1].str());
            std::cout << "Found σres_U2: " << u_res2 << " um" << std::endl;
        }
 
        if (std::regex_search(output, match, u3_pattern)) {
            u_res3 = std::stod(match[1].str());
            std::cout << "Found σres_U3: " << u_res3 << " um" << std::endl;
        }
 
        if (std::regex_search(output, match, v1_pattern)) {
            v_res1 = std::stod(match[1].str());
            std::cout << "Found σres_V1: " << v_res1 << " um" << std::endl;
        }
 
        if (std::regex_search(output, match, v2_pattern)) {
            v_res2 = std::stod(match[1].str());
            std::cout << "Found σres_V2: " << v_res2 << " um" << std::endl;
        }
 
        if (std::regex_search(output, match, v3_pattern)) {
            v_res3 = std::stod(match[1].str());
            std::cout << "Found σres_V3: " << v_res3 << " um" << std::endl;
        }
 
        return {u_res1, u_res2, u_res3, v_res1, v_res2, v_res3};
    }
 
    std::pair<double, double> getTrackErrors() {
        std::ostringstream cmd;
        cmd << "root -l -q -b 'u_v_err.C(\"" << rootFile << "\")'";
 
        std::cout << "Running u_v_err.C with file: " << rootFile << std::endl;
        FILE* pipe = popen(cmd.str().c_str(), "r");
        std::string output;
        if (pipe) {
            char buffer[256];
            while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                output += buffer;
            }
            pclose(pipe);
        }
        std::cout << "u_v_err.C output:" << std::endl << output << std::endl;
 
        std::regex u_pattern(R"(mean_u_err\s*=\s*([\d.]+)\s*um)");
        std::regex v_pattern(R"(mean_v_err\s*=\s*([\d.]+)\s*um)");
        std::smatch match;
        double mean_u_err = 0.0;
        double mean_v_err = 0.0;
 
        if (std::regex_search(output, match, u_pattern)) {
            mean_u_err = std::stod(match[1].str());
            std::cout << "Found mean_u_err: " << mean_u_err << " um" << std::endl;
        } else {
            std::cerr << "Warning: mean_u_err not found in output!" << std::endl;
        }
 
        if (std::regex_search(output, match, v_pattern)) {
            mean_v_err = std::stod(match[1].str());
            std::cout << "Found mean_v_err: " << mean_v_err << " um" << std::endl;
        } else {
            std::cerr << "Warning: mean_v_err not found in output!" << std::endl;
        }
        return {mean_u_err, mean_v_err};
    }

 
    ResolutionParams calculateNewParams(const ResolutionParams& current,
                                       double u_res1, double u_res2, double u_res3,
                                       double v_res1, double v_res2, double v_res3,
                                       double mean_u_err, double mean_v_err) {
        std::cout << "Input σres_U: " << u_res1 << ", " << u_res2 << ", " << u_res3 << std::endl;
        std::cout << "Input σres_V: " << v_res1 << ", " << v_res2 << ", " << v_res3 << std::endl;
        std::cout << "Input σtrack: U=" << mean_u_err << ", V=" << mean_v_err << std::endl;
 
        std::cout << std::fixed << std::setprecision(5);

        // Check whether σres is valid
        if (u_res1 <= 0 || u_res2 <= 0 || u_res3 <= 0 ||
            v_res1 <= 0 || v_res2 <= 0 || v_res3 <= 0) {
            std::cerr << "Invalid σres values, keeping current parameters" << std::endl;
            return current;
        }
 
        if (mean_u_err <= 0 || mean_v_err <= 0) {
            std::cerr << "Invalid σtrack values, keeping current parameters" << std::endl;
            return current;
        }
 
        // 根据公式 σpixel = √(σres² - σtrack²) 计算三种cluster size的σpixel
        double u_diff1 = u_res1 * u_res1 - mean_u_err * mean_u_err;
        double u_diff2 = u_res2 * u_res2 - mean_u_err * mean_u_err;
        double u_diff3 = u_res3 * u_res3 - mean_u_err * mean_u_err;
 
        double v_diff1 = v_res1 * v_res1 - mean_v_err * mean_v_err;
        double v_diff2 = v_res2 * v_res2 - mean_v_err * mean_v_err;
        double v_diff3 = v_res3 * v_res3 - mean_v_err * mean_v_err;
 
        std::cout << "U_diffs: " << u_diff1 << ", " << u_diff2 << ", " << u_diff3 << std::endl;
        std::cout << "V_diffs: " << v_diff1 << ", " << v_diff2 << ", " << v_diff3 << std::endl;
 
        if (u_diff1 <= 0 || u_diff2 <= 0 || u_diff3 <= 0 ||
            v_diff1 <= 0 || v_diff2 <= 0 || v_diff3 <= 0) {
            std::cerr << "Invalid difference (σres² - σtrack² <= 0), keeping current parameters" << std::endl;
            return current;
        }
 
        double u_pixel1 = std::sqrt(u_diff1);
        double u_pixel2 = std::sqrt(u_diff2);
        double u_pixel3 = std::sqrt(u_diff3);
 
        double v_pixel1 = std::sqrt(v_diff1);
        double v_pixel2 = std::sqrt(v_diff2);
        double v_pixel3 = std::sqrt(v_diff3);
 
        std::cout << "Calculated σpixel_U: " << u_pixel1 << ", " << u_pixel2 << ", " << u_pixel3 << std::endl;
        std::cout << "Calculated σpixel_V: " << v_pixel1 << ", " << v_pixel2 << ", " << v_pixel3 << std::endl;
 
        // Limit the minimum value to avoid parameters being too small
        double min_value = 1.0; // Minimum value=1um
        u_pixel1 = std::max(u_pixel1, min_value);
        u_pixel2 = std::max(u_pixel2, min_value);
        u_pixel3 = std::max(u_pixel3, min_value);
        v_pixel1 = std::max(v_pixel1, min_value);
        v_pixel2 = std::max(v_pixel2, min_value);
        v_pixel3 = std::max(v_pixel3, min_value);
 
        ResolutionParams new_params;
        new_params.u_size1 = u_pixel1;
        new_params.u_size2 = u_pixel2;
        new_params.u_sizeOther = u_pixel3;
        new_params.v_size1 = v_pixel1;
        new_params.v_size2 = v_pixel2;
        new_params.v_sizeOther = v_pixel3;
 
        std::cout << "New U params: " << new_params.u_size1 << ", " << new_params.u_size2 << ", " << new_params.u_sizeOther << std::endl;
        std::cout << "New V params: " << new_params.v_size1 << ", " << new_params.v_size2 << ", " << new_params.v_sizeOther << std::endl;
 
        return new_params;
    }
 
    void runTracking(const ResolutionParams& params, int iteration) {
        std::string iter_config = outputDir + "/config_iter_" + std::to_string(iteration) + ".json";
        std::string iter_root = outputDir + "/result_iter_" + std::to_string(iteration) + ".root";
 
        params.saveToJson(iter_config);
 
        std::ostringstream cmd;
        cmd << "./altelActsTrack -cutChiSquared 13.816 "
            << "-daqFiles " << dataFile << " "
            << "-geometryFile " << geometryFile << " "
            << "-resolutionConfig " << iter_config << " "
            << "-rootFile " << iter_root << " "
            << "-targetIds " << targetId << " "
            << "-eventSkip " << eventSkip << " "
            << "-eventMax " << eventMax;
 
        std::cout << "Running tracking command: " << cmd.str() << std::endl;
 
        rootFile = iter_root;
 
        int result = std::system(cmd.str().c_str());
        if (result != 0) {
            throw std::runtime_error("Tracking failed in iteration " + std::to_string(iteration));
        }
    }
 
private:
    void printIterationTable() {
        std::cout << "\n" << std::string(140, '=') << std::endl;
        std::cout << "ITERATION STATISTICS TABLE" << std::endl;
        std::cout << std::string(140, '=') << std::endl;
        
        // Table head
        std::cout << std::left << std::setw(8) << "Iter"
                  << std::setw(12) << "Pixels"
                  << std::setw(12) << "σres_U1" << std::setw(12) << "σres_U2" << std::setw(12) << "σres_U3"
                  << std::setw(12) << "u_err"
                  << std::setw(12) << "σpixel_U1" << std::setw(12) << "σpixel_U2" << std::setw(12) << "σpixel_U3"
                  << std::setw(12) << "σres_V1" << std::setw(12) << "σres_V2" << std::setw(12) << "σres_V3"
                  << std::setw(12) << "v_err"
                  << std::setw(12) << "σpixel_V1" << std::setw(12) << "σpixel_V2" << std::setw(12) << "σpixel_V3"
                  << std::endl;
        
        std::cout << std::string(140, '-') << std::endl;
        
        // Data row
        for (const auto& result : iterationResults) {
            std::cout << std::left << std::setw(8) << result.iteration_num
                      << std::setw(12) << "1,2,>2"
                      << std::setw(12) << std::fixed << std::setprecision(6) << result.u_res1
                      << std::setw(12) << result.u_res2 << std::setw(12) << result.u_res3
                      << std::setw(12) << result.u_err_mean
                      << std::setw(12) << result.u_pixel1 << std::setw(12) << result.u_pixel2 << std::setw(12) << result.u_pixel3
                      << std::setw(12) << result.v_res1 << std::setw(12) << result.v_res2 << std::setw(12) << result.v_res3
                      << std::setw(12) << result.v_err_mean
                      << std::setw(12) << result.v_pixel1 << std::setw(12) << result.v_pixel2 << std::setw(12) << result.v_pixel3
                      << std::endl;
        }
        
        std::cout << std::string(140, '=') << std::endl;
        std::cout << "Note: 1=single pixel, 2=two pixels, >2=more than two pixels" << std::endl;
        std::cout << std::string(140, '=') << std::endl;
    }
 
public:
    void optimize() {
        std::cout << "Starting parameter optimization..." << std::endl;
        std::cout << "Config file: " << configFile << std::endl;
        std::cout << "Root file: " << rootFile << std::endl;
 
        ResolutionParams current_params;
        current_params.loadFromJson(configFile);
 
        int max_iterations = 50;
        int iteration = 0;
        totalIterations = 0;  // Total number of iterations
 
        while (iteration < max_iterations) {
       // while ( true ) {
            totalIterations++;
            std::cout << "\n=== Iteration " << iteration << " ===" << std::endl;
            std::cout << "Current U params: " << current_params.u_size1 << ", "
                     << current_params.u_size2 << ", " << current_params.u_sizeOther << std::endl;
            std::cout << "Current V params: " << current_params.v_size1 << ", "
                     << current_params.v_size2 << ", " << current_params.v_sizeOther << std::endl;
 
            // Step 1: Run tracking
            runTracking(current_params, iteration);
 
            // Step 2: Obtain the σres and σtrack of the three cluster sizes
            auto [u_res1, u_res2, u_res3, v_res1, v_res2, v_res3] = getResidualValues();
            auto [mean_u_err, mean_v_err] = getTrackErrors();
 
            std::cout << "σres_U: " << u_res1 << ", " << u_res2 << ", " << u_res3 << " um" << std::endl;
            std::cout << "σres_V: " << v_res1 << ", " << v_res2 << ", " << v_res3 << " um" << std::endl;
            std::cout << "σtrack: U=" << mean_u_err << " um, V=" << mean_v_err << " um" << std::endl;
 
            // Step 3: Calculate the σpixel value and record the result
            double u_pixel1 = 0.0, u_pixel2 = 0.0, u_pixel3 = 0.0;
            double v_pixel1 = 0.0, v_pixel2 = 0.0, v_pixel3 = 0.0;
            
            //Calculate the σpixel value
            if (u_res1 > 0 && mean_u_err > 0 && (u_res1*u_res1 - mean_u_err*mean_u_err) > 0) {
                u_pixel1 = std::sqrt(u_res1*u_res1 - mean_u_err*mean_u_err);
            }
            if (u_res2 > 0 && mean_u_err > 0 && (u_res2*u_res2 - mean_u_err*mean_u_err) > 0) {
                u_pixel2 = std::sqrt(u_res2*u_res2 - mean_u_err*mean_u_err);
            }
            if (u_res3 > 0 && mean_u_err > 0 && (u_res3*u_res3 - mean_u_err*mean_u_err) > 0) {
                u_pixel3 = std::sqrt(u_res3*u_res3 - mean_u_err*mean_u_err);
            }
            if (v_res1 > 0 && mean_v_err > 0 && (v_res1*v_res1 - mean_v_err*mean_v_err) > 0) {
                v_pixel1 = std::sqrt(v_res1*v_res1 - mean_v_err*mean_v_err);
            }
            if (v_res2 > 0 && mean_v_err > 0 && (v_res2*v_res2 - mean_v_err*mean_v_err) > 0) {
                v_pixel2 = std::sqrt(v_res2*v_res2 - mean_v_err*mean_v_err);
            }
            if (v_res3 > 0 && mean_v_err > 0 && (v_res3*v_res3 - mean_v_err*mean_v_err) > 0) {
                v_pixel3 = std::sqrt(v_res3*v_res3 - mean_v_err*mean_v_err);
            }
            
            // Record the result of this iteration
            IterationResult result;
            result.iteration_num = iteration;
            result.u_res1 = u_res1; result.u_res2 = u_res2; result.u_res3 = u_res3;
            result.v_res1 = v_res1; result.v_res2 = v_res2; result.v_res3 = v_res3;
            result.u_err_mean = mean_u_err; result.v_err_mean = mean_v_err;
            result.u_pixel1 = u_pixel1; result.u_pixel2 = u_pixel2; result.u_pixel3 = u_pixel3;
            result.v_pixel1 = v_pixel1; result.v_pixel2 = v_pixel2; result.v_pixel3 = v_pixel3;
            
            iterationResults.push_back(result);
 
            // Step 4: Calculate new parameters
            ResolutionParams new_params = calculateNewParams(current_params,
                                                             u_res1, u_res2, u_res3,
                                                             v_res1, v_res2, v_res3,
                                                             mean_u_err, mean_v_err);
 
            // Step 5: Check the convergence
            if (current_params.isConverged(new_params, convergenceThreshold)) {
                std::cout << "\nConvergence achieved!" << std::endl;
                current_params = new_params;
                break;
            }
 
            current_params = new_params;
            iteration++;
        }
 
        // Save the final result
        current_params.saveToJson(configFile + "_optimized");
        std::cout << "Total iterations performed: " << totalIterations << std::endl;
        std::cout << "\nOptimization completed!" << std::endl;
        std::cout << "Final U params: " << current_params.u_size1 << ", "
                 << current_params.u_size2 << ", " << current_params.u_sizeOther << std::endl;
        std::cout << "Final V params: " << current_params.v_size1 << ", "
                 << current_params.v_size2 << ", " << current_params.v_sizeOther << std::endl;
                
        // Output statistical table
        printIterationTable();
    }
};



int main(int argc, char* argv[]) {
    // At least 4 parameters are needed：data_file(At least 1)、geometry_file、config_file、output_root
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <data_file1> [data_file2] ... <geometry_file> <config_file> <output_root> [target_id] [max_events] [conv_threshold]" << std::endl;
        std::cerr << "Example: " << argv[0] << " data1.raw data2.raw geo.json config.json output.root 32 100000 0.000001" << std::endl;
        return 1;
    }

    try {
        // Collect all data file paths (the first N parameters, up to non-data file parameters)
        std::vector<std::string> dataPaths;
        int i = 1;  // Starting from the first parameter (argv[0] is the program name)

        // The data file must be followed in sequence as follows：geometry_file(.json)、config_file(.json)、output_root(.root)
        // Collect all data files first (until the first.json file is encountered, which is geometry_file)
        while (i < argc) {
            std::string arg = argv[i];
            // If a.json file is encountered, the data file collection is complete
            if (arg.find(".json") != std::string::npos) {
                break;
            }
            dataPaths.push_back(arg);
            i++;
        }

        if (dataPaths.empty()) {
            std::cerr << "Error: No data files provided!" << std::endl;
            return 1;
        }

        // Analyze subsequent fixed parameters：geometry_file、config_file、output_root
        if (i + 2 >= argc) {  // At least three parameters are required：geometry、config、output
            std::cerr << "Error: Missing required parameters (geometry/config/output)" << std::endl;
            return 1;
        }
        std::string geometryFile = argv[i++];    // geometry_file（first.json）
        std::string configFile = argv[i++];      // config_file（second.json）
        std::string rootFile = argv[i++];        // output_root（.root file）

        // Resolves optional parameters：target_id、max_events、conv_threshold
        int targetId = 32;                       // Default target_id
        if (i < argc) {
            targetId = std::stoi(argv[i++]);
        }

        int maxEvents = 100000;                  // Default maximum number of events
        if (i < argc) {
            maxEvents = std::stoi(argv[i++]);
        }

        double convThreshold = 0.000001;            //Default convergence threshold
        if (i < argc) {
            convThreshold = std::stod(argv[i++]);
        }

        // Print parameter information (for debugging)
        std::cout << "=== Parameters ===" << std::endl;
        std::cout << "Data files (" << dataPaths.size() << "):" << std::endl;
        for (const auto& file : dataPaths) {
            std::cout << "  - " << file << std::endl;
        }
        std::cout << "Geometry file: " << geometryFile << std::endl;
        std::cout << "Config file: " << configFile << std::endl;
        std::cout << "Output root file: " << rootFile << std::endl;
        std::cout << "Target ID: " << targetId << std::endl;
        std::cout << "Max events: " << maxEvents << std::endl;
        std::cout << "Convergence threshold: " << convThreshold << std::endl;

        //Create an optimizer instance (pass the data file vector)
        ParameterOptimizer optimizer(dataPaths, geometryFile, configFile, rootFile, targetId, maxEvents, convThreshold);
        optimizer.optimize();

        // Display the final statistical information
        std::cout << "\n=== OPTIMIZATION SUMMARY ===" << std::endl;
        std::cout << "Total iterations performed: " << optimizer.getTotalIterations() << std::endl;
        std::cout << "Final configuration saved to: " << configFile + "_optimized" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

