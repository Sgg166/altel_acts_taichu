// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <stdexcept>
#include <string>

#include "Acts/EventData/Measurement.hpp"

#include "TelElement.hpp"
#include "myrapidjson.h"
namespace TelActs {

/// Source link class for Alpide pixel hit.
///
/// The source link stores the measuremts, surface
///

class PixelSourceLink;
using TelSourceLink = PixelSourceLink;

class PixelSourceLink {
public:
  PixelSourceLink(const Acts::Surface &surface, Acts::Vector2D values,
                  Acts::BoundMatrix cov)
      : m_values(values), m_cov(cov), m_surface(&surface) {
  }
  /// Must be default_constructible to satisfy SourceLinkConcept.
  PixelSourceLink() = default;
  PixelSourceLink(PixelSourceLink &&) = default;
  PixelSourceLink(const PixelSourceLink &) = default;
  PixelSourceLink &operator=(PixelSourceLink &&) = default;
  PixelSourceLink &operator=(const PixelSourceLink &) = default;

  constexpr const Acts::Surface &referenceSurface() const { return *m_surface; }
  Acts::FittableMeasurement<PixelSourceLink> operator*() const {
    return Acts::Measurement<PixelSourceLink, Acts::BoundIndices,
                             Acts::eBoundLoc0, Acts::eBoundLoc1>{
        m_surface->getSharedPtr(), *this, m_cov.topLeftCorner<2, 2>(),
        m_values[0], m_values[1]};
  }

  // reset the covariance
  // should be done by calibrator?
  void setCovariance(const Acts::BoundMatrix &cov) { m_cov = cov; }

  // get the global position
  Acts::Vector3D globalPosition(const Acts::GeometryContext &gctx) const {
    Acts::Vector3D mom(1, 1, 1);
    Acts::Vector3D global = m_surface->localToGlobal(gctx, m_values, mom);
    return global;
  }

  static std::vector<TelSourceLink> CreateSourceLinks(const JsonValue &js, const std::vector<std::shared_ptr<TelElement>> eles){
    std::vector<TelSourceLink> sourcelinks;
    const auto &layers = js["layers"];
    JsonAllocator jsa;
    JsonValue js_hits(rapidjson::kArrayType);
    for (const auto &layer : layers.GetArray()) {
      size_t id_ext = layer["ext"].GetUint();
      for (const auto &hit : layer["hit"].GetArray()) {
        double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
        double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
        JsonValue js_hit(rapidjson::kObjectType);
        js_hit.AddMember("id", id_ext, jsa);
        js_hit.AddMember("x", x_hit, jsa);
        js_hit.AddMember("y", y_hit, jsa);
        js_hits.PushBack(std::move(js_hit), jsa);
      }
    }

    for(const auto &js_hit : js_hits.GetArray()) {
      size_t id = js_hit["id"].GetUint();
      std::shared_ptr<TelElement> ele;
      for(const auto& anEle: eles){
        if(anEle->id() == id){
          ele = anEle;
          break;
        }
      }
      if(ele){
        Acts::Vector2D loc_hit;
        double x = js_hit["x"].GetDouble();
        double y = js_hit["y"].GetDouble();
        loc_hit << x, y;
          //////////// hit data
        Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
        double resX = 5_um;
        double resY = 5_um;
        cov_hit(0, 0) = resX * resX;
        cov_hit(1, 1) = resY * resY;
        sourcelinks.emplace_back(ele->surface(), loc_hit, cov_hit);
      }
    }
    return sourcelinks;
  }

private:
  Acts::Vector2D m_values;
  Acts::BoundMatrix m_cov;
  // need to store pointers to make the object copyable
  const Acts::Surface *m_surface;
  friend bool operator==(const PixelSourceLink &lhs,
                         const PixelSourceLink &rhs) {
    return lhs.m_values.isApprox(rhs.m_values) and
           lhs.m_cov.isApprox(rhs.m_cov);
  }
};

} // namespace Telescope
