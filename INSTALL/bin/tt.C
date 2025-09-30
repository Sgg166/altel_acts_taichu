// altelTTreeUVErrReader.C
// 使用方法：
// 1. 创建 root 初始化文件 rootlogon.C（绝对路径优先）
// {
//   gROOT->ProcessLine(".include /path/to/altel_acts/source/teldata/event/include/");
//   gROOT->ProcessLine(".include /path/to/altel_acts/source/teldata/root/include/");
//   gROOT->ProcessLine(".L /path/to/altel_acts/source/teldata/root/src/TelEventTTreeReader.cpp");
// }
// 2. 在 rootlogon.C 所在目录运行 ROOT 宏
//       >> root -l 'altelTTreeUVErrReader.C("../path/to/your/output.root")'
 
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <numeric>
#include <algorithm>
 
void altelTTreeUVErrReader(const std::string& rootFilePath, const std::string& outputPrefix = "uve_analysis") {
    std::cout << "分析文件: " << rootFilePath << std::endl;
    
    // 打开 ROOT 文件
    TFile *tfile = new TFile(rootFilePath.c_str(), "READ");
    if (!tfile || !tfile->IsOpen()) {
        std::fprintf(stderr, "无法打开 ROOT 文件\n");
        throw;
    }
    
    // 获取事件树
    TTree *pTree = 0;
    tfile->GetObject("eventTree", pTree);
    if (!pTree) {
        std::fprintf(stderr, "无法获取 eventTree\n");
        throw;
    }
    
    // 设置读取器
    altel::TelEventTTreeReader ttreeReader;
    ttreeReader.setTTree(pTree);
 
    size_t totalNumEvents = ttreeReader.numEvents();
    std::cout << "总事件数: " << totalNumEvents << std::endl;
 
    // 统计容器
    std::vector<double> all_u_err;
    std::vector<double> all_v_err;
    std::map<int, std::vector<double>> det_u_err_map;
    std::map<int, std::vector<double>> det_v_err_map;
    std::map<int, int> det_hit_count;
 
    // 创建直方图
    TH1D* h_u_err = new TH1D("h_u_err", "U方向误差分布;U误差 (mm);计数", 100, -0.1, 0.1);
    TH1D* h_v_err = new TH1D("h_v_err", "V方向误差分布;V误差 (mm);计数", 100, -0.1, 0.1);
    TH2D* h_uerr_verr = new TH2D("h_uerr_verr", "U误差 vs V误差;U误差 (mm);V误差 (mm)", 50, -0.1, 0.1, 50, -0.1, 0.1);
    
    // 按检测器创建直方图
    std::map<int, TH1D*> h_u_err_det;
    std::map<int, TH1D*> h_v_err_det;
 
    // 遍历所有事件
    for (size_t eventNum = 0; eventNum < totalNumEvents; eventNum++) {
        if (eventNum % 1000 == 0) {
            std::cout << "处理事件: " << eventNum << "/" << totalNumEvents << std::endl;
        }
 
        try {
            std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
            if (!telEvent) continue;
 
            // 遍历所有轨迹
            for (auto aTraj : telEvent->trajs()) {
                if (!aTraj) continue;
 
                // 跳过拟合点数不足的轨迹
                if (aTraj->numOriginMeasHit() < 3) {
                    continue;
                }
 
                // 遍历轨迹中的所有命中点
                for (auto &aTrajHit : aTraj->trajHits()) {
                    if (!aTrajHit) continue;
 
                    auto aFitHit = aTrajHit->fitHit();
                    if (!aFitHit) continue;
 
                    // 获取误差数据
                    double u_err = aFitHit->u_err();
                    double v_err = aFitHit->v_err();
                    int detId = aFitHit->detN();
 
                    // 收集数据
                    all_u_err.push_back(u_err);
                    all_v_err.push_back(v_err);
                    det_u_err_map[detId].push_back(u_err);
                    det_v_err_map[detId].push_back(v_err);
                    det_hit_count[detId]++;
 
                    // 填充总体直方图
                    h_u_err->Fill(u_err);
                    h_v_err->Fill(v_err);
                    h_uerr_verr->Fill(u_err, v_err);
 
                    // 为每个检测器创建直方图（如果不存在）
                    if (h_u_err_det.find(detId) == h_u_err_det.end()) {
                        std::string h_name_u = "h_u_err_det_" + std::to_string(detId);
                        std::string h_title_u = "U误差分布 - 检测器" + std::to_string(detId) + ";U误差 (mm);计数";
                        h_u_err_det[detId] = new TH1D(h_name_u.c_str(), h_title_u.c_str(), 100, -0.1, 0.1);
 
                        std::string h_name_v = "h_v_err_det_" + std::to_string(detId);
                        std::string h_title_v = "V误差分布 - 检测器" + std::to_string(detId) + ";V误差 (mm);计数";
                        h_v_err_det[detId] = new TH1D(h_name_v.c_str(), h_title_v.c_str(), 100, -0.1, 0.1);
                    }
 
                    // 填充检测器直方图
                    h_u_err_det[detId]->Fill(u_err);
                    h_v_err_det[detId]->Fill(v_err);
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "处理事件 " << eventNum << " 时出错: " << e.what() << std::endl;
            continue;
        }
    }
 
    // 计算统计函数
    auto calculateStats = [](const std::vector<double>& data) {
        if (data.empty()) {
            return std::make_tuple(0.0, 0.0, 0.0, 0.0, 0);
        }
 
        double sum = std::accumulate(data.begin(), data.end(), 0.0);
        double mean = sum / data.size();
        
        double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
        double std_dev = std::sqrt(sq_sum / data.size() - mean * mean);
        
        auto minmax = std::minmax_element(data.begin(), data.end());
        
        return std::make_tuple(mean, std_dev, *minmax.first, *minmax.second, static_cast<int>(data.size()));
    };
 
    // 输出总体统计信息
    std::cout << "\n=== 总体统计信息 ===" << std::endl;
    auto [u_mean, u_std, u_min, u_max, u_count] = calculateStats(all_u_err);
    auto [v_mean, v_std, v_min, v_max, v_count] = calculateStats(all_v_err);
    
    std::cout << "U误差统计:" << std::endl;
    std::cout << "  数据点数: " << u_count << std::endl;
    std::cout << "  均值: " << u_mean << " mm" << std::endl;
    std::cout << "  标准差: " << u_std << " mm" << std::endl;
    std::cout << "  范围: [" << u_min << ", " << u_max << "] mm" << std::endl;
    
    std::cout << "V误差统计:" << std::endl;
    std::cout << "  数据点数: " << v_count << std::endl;
    std::cout << "  均值: " << v_mean << " mm" << std::endl;
    std::cout << "  标准差: " << v_std << " mm" << std::endl;
    std::cout << "  范围: [" << v_min << ", " << v_max << "] mm" << std::endl;
 
    // 按检测器输出统计信息
    std::cout << "\n=== 按检测器统计 ===" << std::endl;
    for (const auto& [detId, u_errs] : det_u_err_map) {
        const auto& v_errs = det_v_err_map[detId];
        auto [det_u_mean, det_u_std, det_u_min, det_u_max, det_u_count] = calculateStats(u_errs);
        auto [det_v_mean, det_v_std, det_v_min, det_v_max, det_v_count] = calculateStats(v_errs);
        
        std::cout << "检测器 " << detId << ":" << std::endl;
        std::cout << "  命中点数: " << det_hit_count[detId] << std::endl;
        std::cout << "  U误差 - 均值: " << det_u_mean << " mm, 标准差: " << det_u_std << " mm" << std::endl;
        std::cout << "  V误差 - 均值: " << det_v_mean << " mm, 标准差: " << det_v_std << " mm" << std::endl;
    }
 
    // 保存统计结果到文本文件
    std::string stats_filename = outputPrefix + "_statistics.txt";
    std::ofstream stats_file(stats_filename);
    stats_file << "U/V误差统计分析结果\n";
    stats_file << "输入文件: " << rootFilePath << "\n\n";
    
    stats_file << "总体统计:\n";
    stats_file << "U误差: 均值=" << u_mean << " mm, 标准差=" << u_std << " mm, 最小值=" << u_min << " mm, 最大值=" << u_max << " mm, 数据点数=" << u_count << "\n";
    stats_file << "V误差: 均值=" << v_mean << " mm, 标准差=" << v_std << " mm, 最小值=" << v_min << " mm, 最大值=" << v_max << " mm, 数据点数=" << v_count << "\n\n";
    
    stats_file << "按检测器统计:\n";
    for (const auto& [detId, u_errs] : det_u_err_map) {
        const auto& v_errs = det_v_err_map[detId];
        auto [det_u_mean, det_u_std, det_u_min, det_u_max, det_u_count] = calculateStats(u_errs);
        auto [det_v_mean, det_v_std, det_v_min, det_v_max, det_v_count] = calculateStats(v_errs);
        
        stats_file << "检测器" << detId << ": ";
        stats_file << "U(均值=" << det_u_mean << ", 标准差=" << det_u_std << "), ";
        stats_file << "V(均值=" << det_v_mean << ", 标准差=" << det_v_std << "), ";
        stats_file << "数据点=" << det_hit_count[detId] << "\n";
    }
    stats_file.close();
 
    // 创建画布并绘制图形
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
 
    // U误差分布
    TCanvas* c1 = new TCanvas("c1", "U误差分布", 800, 600);
    c1->cd();
    h_u_err->Draw();
    h_u_err->Fit("gaus", "Q");
    h_u_err->SetLineColor(kBlue);
    h_u_err->SetLineWidth(2);
    c1->SaveAs((outputPrefix + "_u_err.png").c_str());
 
    // V误差分布
    TCanvas* c2 = new TCanvas("c2", "V误差分布", 800, 600);
    c2->cd();
    h_v_err->Draw();
    h_v_err->Fit("gaus", "Q");
    h_v_err->SetLineColor(kRed);
    h_v_err->SetLineWidth(2);
    c2->SaveAs((outputPrefix + "_v_err.png").c_str());
 
    // U误差 vs V误差
    TCanvas* c3 = new TCanvas("c3", "U误差 vs V误差", 800, 600);
    c3->cd();
    h_uerr_verr->Draw("COLZ");
    h_uerr_verr->SetStats(false);
    c3->SaveAs((outputPrefix + "_uerr_verr.png").c_str());
 
    // 按检测器绘制U误差分布
    TCanvas* c4 = new TCanvas("c4", "U误差分布（按检测器）", 1200, 800);
    c4->Divide(2, 2);
    int pad = 1;
    for (const auto& [detId, hist] : h_u_err_det) {
        if (pad <= 4) { // 只显示前4个检测器
            c4->cd(pad++);
            hist->Draw();
            hist->Fit("gaus", "Q");
            hist->SetLineColor(pad);
            hist->SetLineWidth(2);
        }
    }
    c4->SaveAs((outputPrefix + "_u_err_by_detector.png").c_str());
 
    // 保存 ROOT 文件
    std::string root_output = outputPrefix + "_output.root";
    TFile* output_file = new TFile(root_output.c_str(), "RECREATE");
    h_u_err->Write();
    h_v_err->Write();
    h_uerr_verr->Write();
    for (auto& [detId, hist] : h_u_err_det) {
        hist->Write();
    }
    for (auto& [detId, hist] : h_v_err_det) {
        hist->Write();
    }
    output_file->Close();
 
    // 清理内存
    tfile->Close();
 
    std::cout << "\n分析完成！" << std::endl;
    std::cout << "输出文件:" << std::endl;
    std::cout << "  - " << stats_filename << " (统计结果)" << std::endl;
    std::cout << "  - " << outputPrefix << "_u_err.png (U误差分布)" << std::endl;
    std::cout << "  - " << outputPrefix << "_v_err.png (V误差分布)" << std::endl;
    std::cout << "  - " << outputPrefix << "_uerr_verr.png (U-V误差相关性)" << std::endl;
    std::cout << "  - " << outputPrefix << "_u_err_by_detector.png (按检测器U误差)" << std::endl;
    std::cout << "  - " << root_output << " (ROOT文件)" << std::endl;
 
    return;
}
