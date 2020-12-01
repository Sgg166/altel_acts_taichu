
#include "FEI4Helper.hh"
#include "CvtEudaqAltelRaw.hh"

std::shared_ptr<altel::TelEvent> altel::createTelEvent(eudaq::EventSPC eudaqEvent){

  eudaq::EventSPC ev_altel;
  eudaq::EventSPC ev_fei4;

  uint32_t runN = eudaqEvent->GetRunN();
  uint32_t eventN = eudaqEvent->GetEventN();
  uint32_t triggerN = eudaqEvent->GetTriggerN();
  uint32_t deviceN = eudaqEvent->GetDeviceN();
  std::shared_ptr<altel::TelEvent> telev(new altel::TelEvent(runN, eventN, deviceN, triggerN));

  if(eudaqEvent->IsFlagPacket()){
    auto subev_col = eudaqEvent->GetSubEvents();
    for(auto& subev: subev_col){
      if(subev->GetDescription() == "AltelRaw"){
        ev_altel = subev;
      }
      if(subev->GetDescription() == "USBPIXI4"){
        ev_fei4 = subev;
      }
    }
    subev_col.clear();
  }
  else if(eudaqEvent->GetDescription() == "AltelRaw"){
    ev_altel = eudaqEvent;
  }
  if(!ev_altel){
    return telev;
  }

  if(ev_fei4){
    auto ev_raw = std::dynamic_pointer_cast<const eudaq::RawEvent>(ev_fei4);
    auto block_n_list = ev_raw->GetBlockNumList();
    if(block_n_list.size()>1){
      throw;
    }
    if(!block_n_list.empty()){
      auto data = ev_raw->GetBlock(block_n_list[0]);

      std::vector<std::pair<uint16_t, uint16_t>> uvs = FEI4Helper::GetMeasRawUVs(data);
      std::vector<altel::TelMeasRaw> feiMeasRaws;
      for(auto & [uraw, vraw] : uvs){
        feiMeasRaws.emplace_back(uraw, vraw, 101, triggerN); //fei4 detN=101
      }
      auto feiMeasHits = altel::TelMeasHit::clustering_UVDCus(feiMeasRaws,
                                                              FEI4Helper::pitchU,
                                                              FEI4Helper::pitchV,
                                                              -FEI4Helper::pitchU*(FEI4Helper::numPixelU-1)*0.5,
                                                              -FEI4Helper::pitchV*(FEI4Helper::numPixelV-1)*0.5);

      // telev->measRaws().insert(telev->measRaws().end(), feiMeasRaws.begin(), feiMeasRaws.end());
      telev->measHits().insert(telev->measHits().end(), feiMeasHits.begin(), feiMeasHits.end());
    }
  }

  auto ev_altelraw = std::dynamic_pointer_cast<const eudaq::RawEvent>(ev_altel);
  size_t nblocks= ev_altelraw->NumBlocks();
  auto block_n_list = ev_altelraw->GetBlockNumList();
  if(!nblocks)
    throw;

  for(const auto& blockNum: block_n_list){
    auto rawblock = ev_altelraw->GetBlock(blockNum);
    uint32_t* p_block = reinterpret_cast<uint32_t*>(rawblock.data());
    uint32_t layerID = *p_block;

    p_block++;
    uint32_t cluster_size =  *p_block;
    p_block++;

    std::vector<altel::TelMeasRaw> alpideMeasRaws;
    for(size_t i = 0; i < cluster_size; i++){
      float cluster_x_mm = *(reinterpret_cast<float*>(p_block));
      p_block++;
      float cluster_y_mm = *(reinterpret_cast<float*>(p_block));
      p_block++;
      uint32_t pixel_size =  *p_block;
      p_block++;

      for(size_t j = 0; j < pixel_size; j++){
        uint32_t pixelXY =  *p_block;
        p_block++;
        uint16_t pixelX = static_cast<uint16_t>(pixelXY);
        uint16_t pixelY = static_cast<uint16_t>(pixelXY>>16);
        alpideMeasRaws.emplace_back(pixelX, pixelY, layerID, triggerN);
      }
    }
    auto alpideMeasHits = altel::TelMeasHit::clustering_UVDCus(alpideMeasRaws,
                                                               0.02924,
                                                               0.02688,
                                                               -0.02924*(1024-1)*0.5,
                                                               -0.02688*(512-1)*0.5);

    // telev->measRaws().insert(telev->measRaws().end(), alpideMeasRaws.begin(), alpideMeasRaws.end());
    telev->measHits().insert(telev->measHits().end(), alpideMeasHits.begin(), alpideMeasHits.end());
  }
  return telev;
}
