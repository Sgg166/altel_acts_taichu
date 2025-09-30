#include "TelEvent.hpp"
#include <iostream>
#include <memory>
 
using namespace altel;
 
std::shared_ptr<TelEvent> createSimpleTestEvent() {
    auto event = std::make_shared<TelEvent>();
 
    for (int detN = 1; detN <= 6; detN++) {
        auto measHit = std::make_shared<TelMeasHit>();
        measHit->DN = detN;
        measHit->PLs[0] = detN * 1.0;
        measHit->PLs[1] = detN * 1.0;
        measHit->MRs.push_back(TelMeasRaw(detN * 100, detN * 100, detN, 0));
        event->measHits().push_back(measHit);
    }
 
    auto goodTraj = std::make_shared<TelTrajectory>();
    for (int detN = 1; detN <= 6; detN++) {
        auto fitHit = std::make_shared<TelFitHit>(
            detN,
            detN * 1.0, detN * 1.0,
            detN * 10.0, detN * 10.0, detN * 10.0,
            detN * 0.1, detN * 0.1, detN * 0.1,
            event->measHits()[detN-1]
        );
 
        auto trajHit = std::make_shared<TelTrajHit>(detN, fitHit);
 
        if (detN <= 5) {
            trajHit->MM = event->measHits()[detN-1];
        }
 
        goodTraj->trajHits().push_back(trajHit);
    }
    event->trajs().push_back(goodTraj);
 
    auto badTraj = std::make_shared<TelTrajectory>();
    for (int detN = 1; detN <= 2; detN++) {
        auto fitHit = std::make_shared<TelFitHit>(
            detN,
            detN * 2.0, detN * 2.0,
            detN * 20.0, detN * 20.0, detN * 20.0,
            detN * 0.2, detN * 0.2, detN * 0.2,
            event->measHits()[detN-1]
        );
 
        auto trajHit = std::make_shared<TelTrajHit>(detN, fitHit, event->measHits()[detN-1]);
        badTraj->trajHits().push_back(trajHit);
    }
    event->trajs().push_back(badTraj);
 
    return event;
}
 
void validateSimpleStatistics(std::shared_ptr<TelEvent> event) {
    std::cout << "=== Simple event verification  ===" << std::endl;
 
    int totalGoodTrajectories = 0;
    int totalAllTrajectories = 0;
    int totalMatchedHits = 0;
    int goodtotalTrajHits = 0;
    int totalFitHits = 0;
    int totalOriginMeasHits = 0;
    int alltotalTrajHits = 0;
 
    std::cout << "测量点总数: " << event->measHits().size() << std::endl;
    std::cout << "轨迹总数: " << event->trajs().size() << std::endl;
 
    for (size_t i = 0; i < event->trajs().size(); i++) {
        const auto& traj = event->trajs()[i];
        if(!traj){
          std::cout << " Warning:" << i << ": is null" << std::endl;
          continue;
        }
 
        totalAllTrajectories++;
        size_t numPoints = traj->numTrajHit();
        std::cout << "\nTrajectory " << i+1 << ": " << numPoints << " numTrajHit";
 
        bool isGoodTrajectory = (numPoints >= 4);
        std::cout << (isGoodTrajectory ? " [good Trajectory]" : " [bad Trajectory]") << std::endl;
 
        if (isGoodTrajectory) {
            totalGoodTrajectories++;
        }
 
        for (size_t j = 0; j < traj->numTrajHit(); j++) {
            auto trajHit = traj->trajHits()[j];
            alltotalTrajHits++; 
 
            if (!trajHit) {
                std::cout << " Warning:" << j << ": is null" << std::endl;
                continue;
            }
 
            if (isGoodTrajectory) {
                goodtotalTrajHits++; 
            }
 
            std::cout << " 点" << j << ": ";
 
            if (trajHit->hasFitHit()) {
                const auto& fitHit = trajHit->fitHit();
                if(fitHit){
                    totalFitHits++;
                    std::cout << " has FitHit";
 
                    if (fitHit->originMeasHit()) {
                        totalOriginMeasHits++;
                        std::cout << "(has OriginMeasHits)";
                    }
                }
            }
 
            if (trajHit->hasMatchedMeasHit()) {
                const auto& matchedMeasHit = trajHit->matchedMeasHit();
                if(matchedMeasHit){
                    totalMatchedHits++;
                    std::cout << " 有MatchedMeasHit";
                }
            }
 
            std::cout << std::endl;
        }
    }
 
    std::cout << "\n=== Statistical results  ===" << std::endl;
    std::cout << "total all Trajectories: " << totalAllTrajectories << std::endl;
    std::cout << "good  Trajectories: " << totalGoodTrajectories << std::endl;
    std::cout << "total TrajHits: " << alltotalTrajHits << std::endl;
    std::cout << "good  total TrajHits: " << goodtotalTrajHits << std::endl;
    std::cout << "total FitHits: " << totalFitHits << std::endl;
    std::cout << "total MatchedHits: " << totalMatchedHits << std::endl;
    std::cout << "total OriginMeasHits: " << totalOriginMeasHits << std::endl;
 
    std::cout << "\n=== Expected results ===" << std::endl;
    std::cout << " total all Trajectories=2 " << std::endl;
    std::cout << " good  Trajectories=1 " <<std::endl;
    std::cout << " total TrajHits=8 (6+2)" << std::endl;
    std::cout << " good  total TrajHits=6 " << std::endl;
    std::cout << " total FitHits=8 (Each TrajHits has the fitHit )" << std::endl;
    std::cout << " total MatchedHits=7 (good tTrajectories has MatchedHits 5+ bad tTrajectories has MatchedHits2)" << std::endl;
    std::cout << " total OriginMeasHits:=8 (Each fitHit has the original OriginMeasHits)" << std::endl;
}
 
int main() {
    auto event = createSimpleTestEvent();
    validateSimpleStatistics(event);
    return 0;
}
