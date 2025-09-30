#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"
#include "eudaq/FileReader.hh"
#include "CvtEudaqAltelRaw.hh"

#include <cstdio>
#include <set>
#include <vector>
#include <string>
#include <memory>
#include <regex>
#include <chrono>

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -eventSkip      <INT>             number of events to skip
  -eventMax       <INT>             max number of events to process
  -daqFiles  <<PATH0> [PATH1]...>   paths to input daq data files (input)
  -jsonFile       <PATH>            path to output json file (output)

Example:
./altelConvert2json -daqFiles eudaqRaw/file.raw -jsonFile output.json -eventMax 10000
)";

int main(int argc, char *argv[]) {
    int64_t eventMaxNum = 0;
    int64_t eventSkipNum = 0;
    std::vector<std::string> rawFilePathCol;
    std::string jsonFilePath;
    int do_verbose = 0;

    // getopt parse
    {
        struct option longopts[] = {
            {"help", no_argument, NULL, 'h'},
            {"verbose", no_argument, NULL, 'v'},
            {"eventSkip", required_argument, NULL, 's'},
            {"eventMax", required_argument, NULL, 'm'},
            {"daqFiles", required_argument, NULL, 'f'},
            {"jsonFile", required_argument, NULL, 'b'},
            {0,0,0,0}
        };

        if(argc == 1){
            std::fprintf(stderr, "%s\n", help_usage.c_str());
            return 1;
        }

        int c, longindex;
        opterr = 1;

        while((c = getopt_long_only(argc, argv, "-", longopts, &longindex)) != -1){
            switch(c){
                case 's': eventSkipNum = std::stoul(optarg); break;
                case 'm': eventMaxNum = std::stoul(optarg); break;
                case 'f': {
                    optind--;
                    while(optind < argc && *argv[optind] != '-') {
                        rawFilePathCol.push_back(argv[optind]);
                        optind++;
                    }
                    break;
                }
                case 'b': jsonFilePath = optarg; break;
                case 'v': do_verbose=1; break;
                case 'h': std::fprintf(stdout,"%s\n", help_usage.c_str()); return 0;
                case 0: break;
                default: return 1;
            }
        }
    }

    if(rawFilePathCol.empty() || jsonFilePath.empty()){
        std::fprintf(stderr,"%s\n", help_usage.c_str());
        return 1;
    }

    std::fprintf(stdout, "%zu daqFiles:\n", rawFilePathCol.size());
    for(auto &f: rawFilePathCol) std::fprintf(stdout,"  %s\n", f.c_str());

    FILE* fd = std::fopen(jsonFilePath.c_str(),"w");
    if(!fd){
        std::fprintf(stderr,"Cannot open output JSON file\n");
        return 1;
    }
    std::fprintf(fd,"[\n");
    bool isFirstEvent = true;

    uint32_t rawFileNum=0;
    size_t readEventNum=0, eventNum=0;
    auto tp_start = std::chrono::system_clock::now();

    eudaq::FileReaderUP reader;
    std::unique_ptr<JsonFileDeserializer> jsreader;
    bool is_eudaq_raw = true;
    if(std::regex_match(rawFilePathCol.front(), std::regex("\\S+\\.json"))) is_eudaq_raw=false;

    while(true){
        if(eventNum >= eventMaxNum && eventMaxNum>0) break;
        std::shared_ptr<altel::TelEvent> fullEvent;

        if(is_eudaq_raw){
            if(!reader){
                if(rawFileNum < rawFilePathCol.size()){
                    std::fprintf(stdout,"processing raw file: %s\n", rawFilePathCol[rawFileNum].c_str());
                    reader = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash("native"), rawFilePathCol[rawFileNum]);
                    rawFileNum++;
                } else break;
            }
            auto eudaqEvent = reader->GetNextEvent();
            if(!eudaqEvent){ reader.reset(); continue; }
            readEventNum++;
            if(readEventNum <= eventSkipNum) continue;
            eventNum++;
            fullEvent = altel::createTelEvent(eudaqEvent);
        } else {
            if(!jsreader){
                if(rawFileNum < rawFilePathCol.size()){
                    std::fprintf(stdout,"processing json file: %s\n", rawFilePathCol[rawFileNum].c_str());
                    jsreader.reset(new JsonFileDeserializer(rawFilePathCol[rawFileNum]));
                    rawFileNum++;
                } else break;
            }
            auto evpack = jsreader->getNextJsonDocument();
            if(evpack.IsNull()){ jsreader.reset(); continue; }
            eventNum++;
            fullEvent = TelActs::createTelEvent(evpack,0,eventNum,0);
        }

        if(!isFirstEvent) std::fprintf(fd,",\n"); else isFirstEvent=false;
        std::fprintf(fd,"{\n\"layers\":[\n");

        uint16_t detN_last = 65535;
        bool firstLayer = true;

        for(auto &aMeasHit : fullEvent->measHits()){
            if(aMeasHit->measRaws().empty()) continue;

            uint16_t detN = aMeasHit->detN();
            bool newLayer = (detN != detN_last);

            static bool firstHitInLayer = true;

            if(newLayer){
                if(!firstLayer) std::fprintf(fd,"\n ]}\n");
                if(!firstLayer) std::fprintf(fd,",\n");
                detN_last = detN;
                firstLayer=false;
                std::fprintf(fd," {\"det\":\"%s\",\"ver\":5,\"tri\":%d,\"cnt\":%d,\"ext\":%d,\"hit\":[",
                             (detN==101)?"fei4":"alpide", fullEvent->clkN(), fullEvent->eveN(), detN);
                firstHitInLayer = true;
            }

            if(!firstHitInLayer) std::fprintf(fd,",");
            firstHitInLayer=false;

            std::fprintf(fd,"{\"pos\":[%f,%f,%d],\"pix\":[", aMeasHit->u(), aMeasHit->v(), aMeasHit->detN());

            bool firstPixel=true;
            for(auto &aMeasRaw: aMeasHit->measRaws()){
                if(!firstPixel) std::fprintf(fd,",");
                firstPixel=false;
                std::fprintf(fd,"[%d,%d,%d]", aMeasRaw.u(), aMeasRaw.v(), aMeasRaw.detN());
            }
            std::fprintf(fd,"]}");
        }

        if(!firstLayer) std::fprintf(fd,"\n ]}\n"); // 结束最后一层
        std::fprintf(fd,"]}\n"); // 事件结束
    }

    auto tp_end = std::chrono::system_clock::now();
    std::chrono::duration<double> dur_diff = tp_end-tp_start;
    double time_s = dur_diff.count();

    std::fprintf(fd,"\n]\n"); // JSON数组结束
    std::fclose(fd);

    std::fprintf(stdout,"processed %zu events in %.2f seconds\n", eventNum, time_s);
    return 0;
}

