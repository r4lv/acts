// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"

// #include "Acts/TrackFitting/GainMatrixSmoother.hpp"
// #include "Acts/TrackFitting/GainMatrixUpdater.hpp"

#include "Acts/TrackFitting/Chi2Fitter.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Plugins/BField/ScalableBField.hpp"
// #include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TrackFittingChi2/TrackFittingChi2Algorithm.hpp"

    namespace {

template <typename track_fitter_t>
struct TrackFitterChi2FunctionImpl {
  track_fitter_t trackFitterChi2;

  TrackFitterChi2FunctionImpl(track_fitter_t&& f)
      : trackFitterChi2(std::move(f)) {}

  ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Result operator()(
      const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Options&
          options) const {
    return trackFitterChi2.fit(sourceLinks, initialParameters, options);
  };
};

// template <typename Fitter>
// struct DirectedFitterFunctionImpl {
//   Fitter fitter;
//   DirectedFitterFunctionImpl(Fitter&& f) : fitter(std::move(f)) {}

//   ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
//       const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
//       const ActsExamples::TrackParameters& initialParameters,
//       const ActsExamples::TrackFittingAlgorithm::TrackFitterOptions& options,
//       const std::vector<const Acts::Surface*>& sSequence) const {
//     return fitter.fit(sourceLinks, initialParameters, options, sSequence);
//   };
// };
}  // namespace

ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Function
ActsExamples::TrackFittingChi2Algorithm::makeTrackFitterChi2Function(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    Options::BFieldVariant magneticField) {
//   using Updater = Acts::GainMatrixUpdater;
//   using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
      [trackingGeometry](auto&& inputField) -> TrackFitterChi2Function {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::EigenStepper<MagneticField>;
        using Navigator = Acts::Navigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using Fitter = Acts::Chi2Fitter<Propagator>;

        // construct all components for the fitter
        MagneticField field(std::move(inputField));
        Stepper stepper(std::move(field));
        Navigator navigator(trackingGeometry);
        navigator.resolvePassive = false;
        navigator.resolveMaterial = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Fitter trackFitterChi2(std::move(propagator));

        // build the fitter functions. owns the fitter object.
        return TrackFitterChi2FunctionImpl<Fitter>(std::move(trackFitterChi2));
      },
      std::move(magneticField));
}
