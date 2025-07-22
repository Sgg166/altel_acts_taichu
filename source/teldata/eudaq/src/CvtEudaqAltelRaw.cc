#include "CvtEudaqAltelRaw.hh"

std::shared_ptr<altel::TelEvent> altel::createTelEvent(eudaq::EventSPC eudaqEvent){

  eudaq::EventSPC ev_altel;

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
    }
    subev_col.clear();
  }
  else if(eudaqEvent->GetDescription() == "AltelRaw"){
    ev_altel = eudaqEvent;
  }
  if(!ev_altel){
    return telev;
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

    std::vector<altel::TelMeasRaw> someMeasRaws;
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
        someMeasRaws.emplace_back(pixelX, pixelY, layerID, triggerN);
      }
    }

    auto someMeasHits = altel::TelMeasHit::clustering_UVDCus(someMeasRaws,
                                                             0.025,
                                                             0.025,
                                                             -0.025*(1024-1)*0.5,
                                                             -0.025*(512-1)*0.5);

    telev->measRaws().insert(telev->measRaws().end(), someMeasRaws.begin(), someMeasRaws.end());
    telev->measHits().insert(telev->measHits().end(), someMeasHits.begin(), someMeasHits.end());
  }
  return telev;
}
