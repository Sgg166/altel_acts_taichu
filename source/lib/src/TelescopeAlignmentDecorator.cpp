// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeAlignmentDecorator.hpp"

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <random>

FW::Telescope::TelescopeAlignmentDecorator::TelescopeAlignmentDecorator(
    const FW::Telescope::TelescopeAlignmentDecorator::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

FW::ProcessCode FW::Telescope::TelescopeAlignmentDecorator::decorate(
    AlgorithmContext& context) {
  // We need to lock the Decorator
  std::lock_guard<std::mutex> alignmentLock(m_alignmentMutex);

  // In which iov batch are we?
  unsigned int iov = context.eventNumber / m_cfg.iovSize;

  // Set the geometry context
  TelescopeDetectorElement::ContextType alignedContext{iov};
  context.geoContext =
      std::make_any<TelescopeDetectorElement::ContextType>(alignedContext);

  return ProcessCode::SUCCESS;
}
