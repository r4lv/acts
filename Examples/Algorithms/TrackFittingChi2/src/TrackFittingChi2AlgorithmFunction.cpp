// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
// #include "Acts/TrackFitting/GainMatrixSmoother.hpp"
// #include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TrackFittingChi2/TrackFittingChi2Algorithm.hpp"
#include "Acts/TrackFitting/Chi2Fitter.hpp"

namespace {

// using Updater = Acts::GainMatrixUpdater;
// using Smoother = Acts::GainMatrixSmoother;
using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using Fitter = Acts::Chi2Fitter<Propagator>;
// using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
// using DirectFitter = Acts::Chi2Fitter<DirectPropagator>;

struct TrackFitterChi2FunctionImpl {
  Fitter trackFitter;

  TrackFitterChi2FunctionImpl(Fitter&& f) : trackFitterChi2(std::move(f)) {}

  ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Result operator()(
      const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Options& options)
      const {
    return trackFitterChi2.fit(sourceLinks, initialParameters, options);
  };
};

ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Function
ActsExamples::TrackFittingChi2Algorithm::makeTrackFitterChi2Function(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField) {
  Stepper stepper(std::move(magneticField));
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitterChi2(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  return TrackFitterChi2FunctionImpl(std::move(trackFitterChi2));
}
