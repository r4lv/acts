// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFitting/Chi2Fitter.hpp"  // RR
// #include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

class TrackFittingChi2Algorithm final : public BareAlgorithm {
 public:
  /// Track fitter function that takes input measurements, initial trackstate
  /// and fitter options and returns some track-fitter-specific result.
  using TrackFitterChi2Options =
      Acts::Chi2FitterOptions<MeasurementCalibrator, Acts::VoidOutlierFinder>;
  using TrackFitterChi2Result =
      Acts::Result<Acts::Chi2FitterResult<IndexSourceLink>>;
  using TrackFitterChi2Function = std::function<TrackFitterChi2Result(
      const std::vector<IndexSourceLink>&, const TrackParameters&,
      const TrackFitterChi2Options&)>;

  /// Fit function that takes the above parameters plus a sorted surface
  /// sequence for the DirectNavigator to follow
//   using DirectedTrackFitterFunction = std::function<TrackFitterResult(
//       const std::vector<IndexSourceLink>&, const TrackParameters&,
//       const TrackFitterOptions&, const std::vector<const Acts::Surface*>&)>;

  /// Create the track fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static TrackFitterChi2Function makeTrackFitterChi2Function(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      Options::BFieldVariant magneticField);

//   static DirectedTrackFitterFunction makeTrackFitterFunction(
//       Options::BFieldVariant magneticField);

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Boolean determining to use DirectNavigator or standard Navigator
    // bool directNavigation; // TODO: remove
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output fitted trajectories collection.
    std::string outputTrajectories;
    /// Type erased fitter function.
    unsigned int nUpdates;
    /// number of update steps
    TrackFitterChi2Function fit;
    /// Type erased direct navigation fitter function
    // DirectedTrackFitterFunction dFit;
    /// Tracking geometry for surface lookup
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  TrackFittingChi2Algorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  /// Helper function to call correct FitterFunction
  TrackFitterChi2Result fitTrack(
      const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const TrackFitterChi2Options& options,
      const std::vector<const Acts::Surface*>& surfSequence) const;

  Config m_cfg;
};

inline ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Result
ActsExamples::TrackFittingChi2Algorithm::fitTrack(
    const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
    const ActsExamples::TrackParameters& initialParameters,
    const Acts::Chi2FitterOptions<MeasurementCalibrator,
                                  Acts::VoidOutlierFinder>& options,
    const std::vector<const Acts::Surface*>& surfSequence) const {
//   if (m_cfg.directNavigation) {
//     return m_cfg.dFit(sourceLinks, initialParameters, options, surfSequence);
//   }

  return m_cfg.fit(sourceLinks, initialParameters, options);
}

}  // namespace ActsExamples
