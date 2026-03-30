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
#include <tuple>
#include <stdexcept>
#include <algorithm>

#include "myrapidjson.h"

class ParameterOptimizer {
public:
    struct DataChunk {
        int eventSkip = 0;
        int eventMax  = 0;
    };

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

    std::vector<DataChunk> dataChunks;

    enum class OptimizeMode {
        OVERALL,
        SPLIT,
        BOTH
    };

    OptimizeMode mode = OptimizeMode::BOTH;

    struct IterationResult {
        int iteration_num = 0;

        double u_res_overall = 0.0;
        double v_res_overall = 0.0;
        double u_pixel_overall = 0.0;
        double v_pixel_overall = 0.0;

        double u_res1 = 0.0, u_res2 = 0.0, u_res3 = 0.0;
        double v_res1 = 0.0, v_res2 = 0.0, v_res3 = 0.0;

        double u_err_mean = 0.0;
        double v_err_mean = 0.0;

        double u_pixel1 = 0.0, u_pixel2 = 0.0, u_pixel3 = 0.0;
        double v_pixel1 = 0.0, v_pixel2 = 0.0, v_pixel3 = 0.0;
    };

    std::vector<IterationResult> iterationResults;

    struct ResolutionParams {
        double u_overall;
        double v_overall;

        double u_size1, u_size2, u_sizeOther;
        double v_size1, v_size2, v_sizeOther;

        ResolutionParams()
            : u_overall(6.0), v_overall(6.0),
              u_size1(6.0), u_size2(6.0), u_sizeOther(6.0),
              v_size1(6.0), v_size2(6.0), v_sizeOther(6.0) {}

        void loadFromJson(const std::string& filename) {
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Cannot open config file: " << filename << std::endl;
                return;
            }

            std::string content((std::istreambuf_iterator<char>(file)),
                                std::istreambuf_iterator<char>());
            file.close();

            JsonDocument doc;
            if (doc.Parse(content.c_str()).HasParseError()) {
                std::cerr << "JSON parse error in file: " << filename << std::endl;
                return;
            }

            std::cout << "Loading config from: " << filename << std::endl;

            if (doc.HasMember("resolutions")) {
                const auto& res = doc["resolutions"];

                if (res.HasMember("u_overall")) {
                    if (res["u_overall"].IsString())
                        u_overall = std::stod(res["u_overall"].GetString());
                    else if (res["u_overall"].IsNumber())
                        u_overall = res["u_overall"].GetDouble();
                }

                if (res.HasMember("v_overall")) {
                    if (res["v_overall"].IsString())
                        v_overall = std::stod(res["v_overall"].GetString());
                    else if (res["v_overall"].IsNumber())
                        v_overall = res["v_overall"].GetDouble();
                }

                if (res.HasMember("u_direction")) {
                    const auto& u = res["u_direction"];
                    if (u.HasMember("size_1")) {
                        if (u["size_1"].IsString()) u_size1 = std::stod(u["size_1"].GetString());
                        else if (u["size_1"].IsNumber()) u_size1 = u["size_1"].GetDouble();
                    }
                    if (u.HasMember("size_2")) {
                        if (u["size_2"].IsString()) u_size2 = std::stod(u["size_2"].GetString());
                        else if (u["size_2"].IsNumber()) u_size2 = u["size_2"].GetDouble();
                    }
                    if (u.HasMember("size_other")) {
                        if (u["size_other"].IsString()) u_sizeOther = std::stod(u["size_other"].GetString());
                        else if (u["size_other"].IsNumber()) u_sizeOther = u["size_other"].GetDouble();
                    }
                }

                if (res.HasMember("v_direction")) {
                    const auto& v = res["v_direction"];
                    if (v.HasMember("size_1")) {
                        if (v["size_1"].IsString()) v_size1 = std::stod(v["size_1"].GetString());
                        else if (v["size_1"].IsNumber()) v_size1 = v["size_1"].GetDouble();
                    }
                    if (v.HasMember("size_2")) {
                        if (v["size_2"].IsString()) v_size2 = std::stod(v["size_2"].GetString());
                        else if (v["size_2"].IsNumber()) v_size2 = v["size_2"].GetDouble();
                    }
                    if (v.HasMember("size_other")) {
                        if (v["size_other"].IsString()) v_sizeOther = std::stod(v["size_other"].GetString());
                        else if (v["size_other"].IsNumber()) v_sizeOther = v["size_other"].GetDouble();
                    }
                }
            }

            std::cout << "u_overall = " << u_overall << "\n";
            std::cout << "v_overall = " << v_overall << "\n";
            std::cout << "u_size1 = " << u_size1 << ", u_size2 = " << u_size2
                      << ", u_sizeOther = " << u_sizeOther << "\n";
            std::cout << "v_size1 = " << v_size1 << ", v_size2 = " << v_size2
                      << ", v_sizeOther = " << v_sizeOther << "\n";
        }

        void saveToJson(const std::string& filename) const {
            JsonDocument doc;
            doc.SetObject();
            JsonAllocator& allocator = doc.GetAllocator();

            JsonValue resolutions(rapidjson::kObjectType);
            JsonValue u_direction(rapidjson::kObjectType);
            JsonValue v_direction(rapidjson::kObjectType);

            resolutions.AddMember("u_overall",
                                  JsonValue().SetString(std::to_string(u_overall).c_str(), allocator),
                                  allocator);
            resolutions.AddMember("v_overall",
                                  JsonValue().SetString(std::to_string(v_overall).c_str(), allocator),
                                  allocator);

            u_direction.AddMember("size_1",
                                  JsonValue().SetString(std::to_string(u_size1).c_str(), allocator),
                                  allocator);
            u_direction.AddMember("size_2",
                                  JsonValue().SetString(std::to_string(u_size2).c_str(), allocator),
                                  allocator);
            u_direction.AddMember("size_other",
                                  JsonValue().SetString(std::to_string(u_sizeOther).c_str(), allocator),
                                  allocator);

            v_direction.AddMember("size_1",
                                  JsonValue().SetString(std::to_string(v_size1).c_str(), allocator),
                                  allocator);
            v_direction.AddMember("size_2",
                                  JsonValue().SetString(std::to_string(v_size2).c_str(), allocator),
                                  allocator);
            v_direction.AddMember("size_other",
                                  JsonValue().SetString(std::to_string(v_sizeOther).c_str(), allocator),
                                  allocator);

            resolutions.AddMember("u_direction", u_direction, allocator);
            resolutions.AddMember("v_direction", v_direction, allocator);
            doc.AddMember("resolutions", resolutions, allocator);

            std::ofstream file(filename);
            rapidjson::StringBuffer buffer;
            rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
            doc.Accept(writer);
            file << buffer.GetString();
        }

        bool isConverged(const ResolutionParams& other, double threshold, OptimizeMode mode) const {
            bool overall_ok = true;
            bool split_ok = true;

            if (mode == OptimizeMode::OVERALL || mode == OptimizeMode::BOTH) {
                overall_ok =
                    std::abs(u_overall - other.u_overall) < threshold &&
                    std::abs(v_overall - other.v_overall) < threshold;
            }

            if (mode == OptimizeMode::SPLIT || mode == OptimizeMode::BOTH) {
                split_ok =
                    std::abs(u_size1 - other.u_size1) < threshold &&
                    std::abs(u_size2 - other.u_size2) < threshold &&
                    std::abs(u_sizeOther - other.u_sizeOther) < threshold &&
                    std::abs(v_size1 - other.v_size1) < threshold &&
                    std::abs(v_size2 - other.v_size2) < threshold &&
                    std::abs(v_sizeOther - other.v_sizeOther) < threshold;
            }

            return overall_ok && split_ok;
        }
    };

    struct OverallResiduals {
        double u_res = 0.0;
        double v_res = 0.0;
    };

    struct SplitResiduals {
        double u_res1 = 0.0, u_res2 = 0.0, u_res3 = 0.0;
        double v_res1 = 0.0, v_res2 = 0.0, v_res3 = 0.0;
    };

public:
    ParameterOptimizer(const std::vector<std::string>& dataPaths,
                       const std::string& geomPath,
                       const std::string& configPath,
                       const std::string& outputPath,
                       int target,
                       int maxEvents = 100000,
                       double convThreshold = 0.000001,
                       const std::string& modeStr = "both",
                       const std::vector<DataChunk>& chunks = {})
        : geometryFile(geomPath), configFile(configPath),
          rootFile(outputPath), targetId(target),
          eventSkip(0), eventMax(maxEvents),
          convergenceThreshold(convThreshold),
          dataChunks(chunks) {

        for (size_t i = 0; i < dataPaths.size(); i++) {
            if (i > 0) dataFile += " ";
            dataFile += dataPaths[i];
        }

        outputDir = std::filesystem::path(configPath).parent_path().string() + "/optimization_results";
        std::filesystem::create_directories(outputDir);

        if (modeStr == "overall") mode = OptimizeMode::OVERALL;
        else if (modeStr == "split") mode = OptimizeMode::SPLIT;
        else mode = OptimizeMode::BOTH;
    }

    int getTotalIterations() const {
        return totalIterations;
    }

private:
    static bool extractValue(const std::string& text, const std::regex& pattern, double& value) {
        std::smatch match;
        if (std::regex_search(text, match, pattern)) {
            value = std::stod(match[1].str());
            return true;
        }
        return false;
    }

    static double calcPixelSigma(double sigmaRes, double sigmaTrack, double minValue = 1.0) {
        if (sigmaRes <= 0 || sigmaTrack <= 0) return 0.0;
        double diff = sigmaRes * sigmaRes - sigmaTrack * sigmaTrack;
        if (diff <= 0) return 0.0;
        return std::max(std::sqrt(diff), minValue);
    }

    std::string runCommandAndCapture(const std::string& cmd) {
        std::cout << "Running command: " << cmd << std::endl;

        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            throw std::runtime_error("Failed to execute command: " + cmd);
        }

        std::string output;
        char buffer[512];
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            output += buffer;
        }

        int ret = pclose(pipe);
        if (ret == -1) {
            throw std::runtime_error("Failed to close pipe for command: " + cmd);
        }

        return output;
    }

    std::string quotePath(const std::string& path) const {
        return "\"" + path + "\"";
    }

    void runOneTrackingJob(const std::string& configPath,
                           const std::string& outRoot,
                           int skip,
                           int maxEvt) {
        std::ostringstream cmd;
        cmd << "./altelActsTrack -cutChiSquared 13.816 "
            << "-daqFiles " << dataFile << " "
            << "-geometryFile " << quotePath(geometryFile) << " "
            << "-resolutionConfig " << quotePath(configPath) << " "
            << "-rootFile " << quotePath(outRoot) << " "
            << "-targetIds " << targetId << " "
            << "-eventSkip " << skip << " "
            << "-eventMax " << maxEvt;

        std::cout << "Running tracking command: " << cmd.str() << std::endl;
        int result = std::system(cmd.str().c_str());
        if (result != 0) {
            throw std::runtime_error("Tracking chunk failed: skip=" + std::to_string(skip) +
                                     ", max=" + std::to_string(maxEvt));
        }
    }

    void mergeRootFiles(const std::vector<std::string>& inputRoots, const std::string& mergedRoot) {
        if (inputRoots.empty()) {
            throw std::runtime_error("No ROOT files to merge.");
        }

        if (inputRoots.size() == 1) {
            std::filesystem::copy_file(inputRoots[0], mergedRoot,
                                       std::filesystem::copy_options::overwrite_existing);
            return;
        }

        std::ostringstream cmd;
        cmd << "hadd -f " << quotePath(mergedRoot) << " ";
        for (const auto& f : inputRoots) {
            cmd << quotePath(f) << " ";
        }

        std::cout << "Running merge command: " << cmd.str() << std::endl;
        int result = std::system(cmd.str().c_str());
        if (result != 0) {
            throw std::runtime_error("hadd failed for merged file: " + mergedRoot);
        }
    }

    OverallResiduals getOverallResidualValues() {
        OverallResiduals result;
        std::ostringstream cmd;
        cmd << "root -l -q -b 'match.C(\"" << rootFile << "\")'";

        std::string output = runCommandAndCapture(cmd.str());
        std::cout << "match.C output:\n" << output << std::endl;

        std::regex u_pattern(R"(SigmaU\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex v_pattern(R"(SigmaV\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex u_pattern_alt(R"(SigmaU1\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex v_pattern_alt(R"(SigmaV1\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");

        if (!extractValue(output, u_pattern, result.u_res)) {
            extractValue(output, u_pattern_alt, result.u_res);
        }
        if (!extractValue(output, v_pattern, result.v_res)) {
            extractValue(output, v_pattern_alt, result.v_res);
        }

        std::cout << "Found overall σres_U = " << result.u_res << " um\n";
        std::cout << "Found overall σres_V = " << result.v_res << " um\n";

        return result;
    }

    SplitResiduals getSplitResidualValues() {
        SplitResiduals result;
        std::ostringstream cmd;
        cmd << "root -l -q -b 'Match.C(\"" << rootFile << "\")'";

        std::string output = runCommandAndCapture(cmd.str());
        std::cout << "Match.C output:\n" << output << std::endl;

        std::regex u1_pattern(R"(SigmaU1\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex u2_pattern(R"(SigmaU2\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex u3_pattern(R"(SigmaU3\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");

        std::regex v1_pattern(R"(SigmaV1\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex v2_pattern(R"(SigmaV2\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex v3_pattern(R"(SigmaV3\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");

        extractValue(output, u1_pattern, result.u_res1);
        extractValue(output, u2_pattern, result.u_res2);
        extractValue(output, u3_pattern, result.u_res3);

        extractValue(output, v1_pattern, result.v_res1);
        extractValue(output, v2_pattern, result.v_res2);
        extractValue(output, v3_pattern, result.v_res3);

        std::cout << "Found split σres_U = "
                  << result.u_res1 << ", " << result.u_res2 << ", " << result.u_res3 << " um\n";
        std::cout << "Found split σres_V = "
                  << result.v_res1 << ", " << result.v_res2 << ", " << result.v_res3 << " um\n";

        return result;
    }

    std::pair<double, double> getTrackErrors() {
        std::ostringstream cmd;
        cmd << "root -l -q -b 'u_v_err.C(\"" << rootFile << "\")'";

        std::string output = runCommandAndCapture(cmd.str());
        std::cout << "u_v_err.C output:\n" << output << std::endl;

        std::regex u_pattern(R"(mean_u_err\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");
        std::regex v_pattern(R"(mean_v_err\s*=\s*([-+]?[\d]*\.?[\d]+)\s*um)");

        double mean_u_err = 0.0;
        double mean_v_err = 0.0;

        if (!extractValue(output, u_pattern, mean_u_err)) {
            std::cerr << "Warning: mean_u_err not found in output!" << std::endl;
        }
        if (!extractValue(output, v_pattern, mean_v_err)) {
            std::cerr << "Warning: mean_v_err not found in output!" << std::endl;
        }

        std::cout << "Found mean_u_err = " << mean_u_err << " um\n";
        std::cout << "Found mean_v_err = " << mean_v_err << " um\n";

        return {mean_u_err, mean_v_err};
    }

    ResolutionParams calculateNewParams(const ResolutionParams& current,
                                        const OverallResiduals& overallRes,
                                        const SplitResiduals& splitRes,
                                        double mean_u_err,
                                        double mean_v_err) {
        ResolutionParams new_params = current;

        std::cout << std::fixed << std::setprecision(6);

        if (mode == OptimizeMode::OVERALL || mode == OptimizeMode::BOTH) {
            double u_pixel = calcPixelSigma(overallRes.u_res, mean_u_err);
            double v_pixel = calcPixelSigma(overallRes.v_res, mean_v_err);

            if (u_pixel > 0) new_params.u_overall = u_pixel;
            if (v_pixel > 0) new_params.v_overall = v_pixel;

            std::cout << "New overall params: U=" << new_params.u_overall
                      << ", V=" << new_params.v_overall << std::endl;
        }

        if (mode == OptimizeMode::SPLIT || mode == OptimizeMode::BOTH) {
            double u_pixel1 = calcPixelSigma(splitRes.u_res1, mean_u_err);
            double u_pixel2 = calcPixelSigma(splitRes.u_res2, mean_u_err);
            double u_pixel3 = calcPixelSigma(splitRes.u_res3, mean_u_err);

            double v_pixel1 = calcPixelSigma(splitRes.v_res1, mean_v_err);
            double v_pixel2 = calcPixelSigma(splitRes.v_res2, mean_v_err);
            double v_pixel3 = calcPixelSigma(splitRes.v_res3, mean_v_err);

            if (u_pixel1 > 0) new_params.u_size1 = u_pixel1;
            if (u_pixel2 > 0) new_params.u_size2 = u_pixel2;
            if (u_pixel3 > 0) new_params.u_sizeOther = u_pixel3;

            if (v_pixel1 > 0) new_params.v_size1 = v_pixel1;
            if (v_pixel2 > 0) new_params.v_size2 = v_pixel2;
            if (v_pixel3 > 0) new_params.v_sizeOther = v_pixel3;

            std::cout << "New split U params: "
                      << new_params.u_size1 << ", "
                      << new_params.u_size2 << ", "
                      << new_params.u_sizeOther << std::endl;

            std::cout << "New split V params: "
                      << new_params.v_size1 << ", "
                      << new_params.v_size2 << ", "
                      << new_params.v_sizeOther << std::endl;
        }

        return new_params;
    }

public:
    void runTracking(const ResolutionParams& params, int iteration) {
        std::string iter_config = outputDir + "/config_iter_" + std::to_string(iteration) + ".json";
        params.saveToJson(iter_config);

        if (dataChunks.empty()) {
            std::string iter_root = outputDir + "/result_iter_" + std::to_string(iteration) + ".root";
            runOneTrackingJob(iter_config, iter_root, eventSkip, eventMax);
            rootFile = iter_root;
            return;
        }

        std::vector<std::string> chunkRoots;
        for (size_t i = 0; i < dataChunks.size(); ++i) {
            std::string part_root = outputDir + "/result_iter_" + std::to_string(iteration)
                                  + "_part_" + std::to_string(i) + ".root";

            runOneTrackingJob(iter_config,
                              part_root,
                              dataChunks[i].eventSkip,
                              dataChunks[i].eventMax);

            chunkRoots.push_back(part_root);
        }

        std::string merged_root = outputDir + "/result_iter_" + std::to_string(iteration) + "_merged.root";
        mergeRootFiles(chunkRoots, merged_root);
        rootFile = merged_root;
    }

private:
    void printIterationTable() {
        std::cout << "\n" << std::string(220, '=') << std::endl;
        std::cout << "ITERATION STATISTICS TABLE" << std::endl;
        std::cout << std::string(220, '=') << std::endl;

        std::cout << std::left
                  << std::setw(6)  << "Iter"
                  << std::setw(12) << "Ures"
                  << std::setw(12) << "Ures1"
                  << std::setw(12) << "Ures2"
                  << std::setw(12) << "Ures3"
                  << std::setw(12) << "Uerr"
                  << std::setw(12) << "Upix"
                  << std::setw(12) << "Upix1"
                  << std::setw(12) << "Upix2"
                  << std::setw(12) << "Upix3"
                  << std::setw(12) << "Vres"
                  << std::setw(12) << "Vres1"
                  << std::setw(12) << "Vres2"
                  << std::setw(12) << "Vres3"
                  << std::setw(12) << "Verr"
                  << std::setw(12) << "Vpix"
                  << std::setw(12) << "Vpix1"
                  << std::setw(12) << "Vpix2"
                  << std::setw(12) << "Vpix3"
                  << std::endl;

        std::cout << std::string(220, '-') << std::endl;

        for (const auto& r : iterationResults) {
            std::cout << std::left
                      << std::setw(6)  << r.iteration_num
                      << std::setw(12) << std::fixed << std::setprecision(6) << r.u_res_overall
                      << std::setw(12) << r.u_res1
                      << std::setw(12) << r.u_res2
                      << std::setw(12) << r.u_res3
                      << std::setw(12) << r.u_err_mean
                      << std::setw(12) << r.u_pixel_overall
                      << std::setw(12) << r.u_pixel1
                      << std::setw(12) << r.u_pixel2
                      << std::setw(12) << r.u_pixel3
                      << std::setw(12) << r.v_res_overall
                      << std::setw(12) << r.v_res1
                      << std::setw(12) << r.v_res2
                      << std::setw(12) << r.v_res3
                      << std::setw(12) << r.v_err_mean
                      << std::setw(12) << r.v_pixel_overall
                      << std::setw(12) << r.v_pixel1
                      << std::setw(12) << r.v_pixel2
                      << std::setw(12) << r.v_pixel3
                      << std::endl;
        }

        std::cout << std::string(220, '=') << std::endl;
        std::cout << "Overall: U/V total resolution. Split: 1=single pixel, 2=double pixel, 3=multi pixel(>2)." << std::endl;
        std::cout << std::string(220, '=') << std::endl;
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
        totalIterations = 0;

        while (iteration < max_iterations) {
            totalIterations++;
            std::cout << "\n=== Iteration " << iteration << " ===" << std::endl;

            std::cout << "Current overall params: U=" << current_params.u_overall
                      << ", V=" << current_params.v_overall << std::endl;
            std::cout << "Current split U params: "
                      << current_params.u_size1 << ", "
                      << current_params.u_size2 << ", "
                      << current_params.u_sizeOther << std::endl;
            std::cout << "Current split V params: "
                      << current_params.v_size1 << ", "
                      << current_params.v_size2 << ", "
                      << current_params.v_sizeOther << std::endl;

            runTracking(current_params, iteration);

            OverallResiduals overallRes;
            SplitResiduals splitRes;

            if (mode == OptimizeMode::OVERALL || mode == OptimizeMode::BOTH) {
                overallRes = getOverallResidualValues();
            }

            if (mode == OptimizeMode::SPLIT || mode == OptimizeMode::BOTH) {
                splitRes = getSplitResidualValues();
            }

            auto [mean_u_err, mean_v_err] = getTrackErrors();

            IterationResult result;
            result.iteration_num = iteration;
            result.u_err_mean = mean_u_err;
            result.v_err_mean = mean_v_err;

            if (mode == OptimizeMode::OVERALL || mode == OptimizeMode::BOTH) {
                result.u_res_overall = overallRes.u_res;
                result.v_res_overall = overallRes.v_res;
                result.u_pixel_overall = calcPixelSigma(overallRes.u_res, mean_u_err);
                result.v_pixel_overall = calcPixelSigma(overallRes.v_res, mean_v_err);
            }

            if (mode == OptimizeMode::SPLIT || mode == OptimizeMode::BOTH) {
                result.u_res1 = splitRes.u_res1;
                result.u_res2 = splitRes.u_res2;
                result.u_res3 = splitRes.u_res3;
                result.v_res1 = splitRes.v_res1;
                result.v_res2 = splitRes.v_res2;
                result.v_res3 = splitRes.v_res3;

                result.u_pixel1 = calcPixelSigma(splitRes.u_res1, mean_u_err);
                result.u_pixel2 = calcPixelSigma(splitRes.u_res2, mean_u_err);
                result.u_pixel3 = calcPixelSigma(splitRes.u_res3, mean_u_err);
                result.v_pixel1 = calcPixelSigma(splitRes.v_res1, mean_v_err);
                result.v_pixel2 = calcPixelSigma(splitRes.v_res2, mean_v_err);
                result.v_pixel3 = calcPixelSigma(splitRes.v_res3, mean_v_err);
            }

            iterationResults.push_back(result);

            ResolutionParams new_params = calculateNewParams(
                current_params, overallRes, splitRes, mean_u_err, mean_v_err
            );

            if (current_params.isConverged(new_params, convergenceThreshold, mode)) {
                std::cout << "\nConvergence achieved!" << std::endl;
                current_params = new_params;
                break;
            }

            current_params = new_params;
            iteration++;
        }

        current_params.saveToJson(configFile + "_optimized");

        std::cout << "Total iterations performed: " << totalIterations << std::endl;
        std::cout << "\nOptimization completed!" << std::endl;
        std::cout << "Final overall params: U=" << current_params.u_overall
                  << ", V=" << current_params.v_overall << std::endl;
        std::cout << "Final split U params: "
                  << current_params.u_size1 << ", "
                  << current_params.u_size2 << ", "
                  << current_params.u_sizeOther << std::endl;
        std::cout << "Final split V params: "
                  << current_params.v_size1 << ", "
                  << current_params.v_size2 << ", "
                  << current_params.v_sizeOther << std::endl;

        printIterationTable();
    }
};

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static std::vector<ParameterOptimizer::DataChunk> loadChunksFromFile(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open chunk file: " + filename);
    }

    std::vector<ParameterOptimizer::DataChunk> chunks;
    std::string line;
    int lineNo = 0;

    while (std::getline(fin, line)) {
        ++lineNo;
        line = trim(line);

        if (line.empty()) continue;
        if (line[0] == '#') continue;

        ParameterOptimizer::DataChunk c;
        size_t pos = line.find(':');

        if (pos != std::string::npos) {
            c.eventSkip = std::stoi(trim(line.substr(0, pos)));
            c.eventMax  = std::stoi(trim(line.substr(pos + 1)));
        } else {
            std::istringstream iss(line);
            if (!(iss >> c.eventSkip >> c.eventMax)) {
                throw std::runtime_error("Invalid chunk format in " + filename +
                                         " at line " + std::to_string(lineNo) +
                                         ": " + line);
            }
        }

        chunks.push_back(c);
    }

    if (chunks.empty()) {
        throw std::runtime_error("No valid chunks found in file: " + filename);
    }

    return chunks;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <data_file1> [data_file2] ... <geometry_file.json> <config_file.json> <output_root.root> "
                  << "[target_id] [max_events] [conv_threshold] [mode] "
                  << "[--chunk skip:max ...] [--chunkFile chunks.txt]\n";
        std::cerr << "mode = overall | split | both\n";
        std::cerr << "Examples:\n";
        std::cerr << argv[0]
                  << " data1.raw data2.raw geo.json config.json output.root 32 100000 0.000001 both "
                  << "--chunk 0:100000 --chunk 200000:200000\n";
        std::cerr << argv[0]
                  << " data1.raw data2.raw geo.json config.json output.root 32 100000 0.000001 both "
                  << "--chunkFile chunks.txt\n";
        return 1;
    }

    try {
        std::vector<std::string> dataPaths;
        int i = 1;

        while (i < argc) {
            std::string arg = argv[i];
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

        if (i + 2 >= argc) {
            std::cerr << "Error: Missing required parameters (geometry/config/output)" << std::endl;
            return 1;
        }

        std::string geometryFile = argv[i++];
        std::string configFile   = argv[i++];
        std::string rootFile     = argv[i++];

        int targetId = 32;
        if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) targetId = std::stoi(argv[i++]);

        int maxEvents = 100000;
        if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) maxEvents = std::stoi(argv[i++]);

        double convThreshold = 0.000001;
        if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) convThreshold = std::stod(argv[i++]);

        std::string mode = "both";
        if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) mode = argv[i++];

        std::vector<ParameterOptimizer::DataChunk> chunks;

        while (i < argc) {
            std::string arg = argv[i++];

            if (arg == "--chunk") {
                if (i >= argc) {
                    throw std::runtime_error("Missing value after --chunk");
                }

                std::string spec = argv[i++];
                auto pos = spec.find(':');
                if (pos == std::string::npos) {
                    throw std::runtime_error("Invalid chunk format: " + spec + ", expected skip:max");
                }

                ParameterOptimizer::DataChunk c;
                c.eventSkip = std::stoi(spec.substr(0, pos));
                c.eventMax  = std::stoi(spec.substr(pos + 1));
                chunks.push_back(c);
            }
            else if (arg == "--chunkFile") {
                if (i >= argc) {
                    throw std::runtime_error("Missing filename after --chunkFile");
                }
                std::string chunkFile = argv[i++];
                auto fileChunks = loadChunksFromFile(chunkFile);
                chunks.insert(chunks.end(), fileChunks.begin(), fileChunks.end());
            }
            else {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }

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
        std::cout << "Mode: " << mode << std::endl;

        if (!chunks.empty()) {
            std::cout << "Chunks:" << std::endl;
            for (size_t k = 0; k < chunks.size(); ++k) {
                std::cout << "  chunk " << k
                          << " -> eventSkip=" << chunks[k].eventSkip
                          << ", eventMax=" << chunks[k].eventMax << std::endl;
            }
        }

        ParameterOptimizer optimizer(
            dataPaths,
            geometryFile,
            configFile,
            rootFile,
            targetId,
            maxEvents,
            convThreshold,
            mode,
            chunks
        );

        optimizer.optimize();

        std::cout << "\n=== OPTIMIZATION SUMMARY ===" << std::endl;
        std::cout << "Total iterations performed: " << optimizer.getTotalIterations() << std::endl;
        std::cout << "Final configuration saved to: " << configFile + "_optimized" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
