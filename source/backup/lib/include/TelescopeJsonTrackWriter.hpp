// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <fstream>
#include <mutex>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "PixelMultiTrajectory.hpp"

#include "JsonGenerator.hpp"

namespace Telescope {

/// @class ObjTelescopeTrackWriter
///
/// Write out the tracks reconstructed using Combinatorial Kalman Filter
/// Writes one file per event with form:
/// One Thread per write call and hence thread safe
class TelescopeJsonTrackWriter
    : public ActsExamples::WriterT<std::vector<PixelMultiTrajectory>> {
public:
  struct Config {
    std::string inputTrajectories; ///< input (fitted) trajectories collection
    std::map<size_t, std::shared_ptr<const Acts::Surface>> trackSurfaces;
    std::string outputDir;
    std::string outputFileName;
  };

  /// Constructor with arguments
  ///
  /// @param cfg configuration struct
  /// @param level Output logging level
  TelescopeJsonTrackWriter(const Config &cfg,
                           Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~TelescopeJsonTrackWriter() override = default;

  /// End-of-run hook
  ActsExamples::ProcessCode endRun() final override;

private:
  Config m_cfg; ///!< Internal configuration represenation
  std::map<std::shared_ptr<const Acts::Surface>, size_t> m_surface_id_map;

  std::unique_ptr<Telescope::JsonAllocator> m_jsa;
  char m_jsbuffer[UINT16_MAX];
  std::FILE *m_jsfp;
  std::unique_ptr<rapidjson::FileWriteStream> m_jsos;
  std::unique_ptr<rapidjson::Writer<rapidjson::FileWriteStream>> m_jsw;
  std::mutex m_mtx;

protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ActsExamples::ProcessCode writeT(
      const ActsExamples::AlgorithmContext &context,
      const std::vector<PixelMultiTrajectory> &trackCollection) final override;
};
} // namespace Telescope
