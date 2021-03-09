// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

/// Combined options for the chi square fitter.
///
/// @tparam calibrator_t Source link type, should be semiregular.
/// @tparam outlier_finder_t Outlier finder type, shoule be semiregular.
template <typename calibrator_t, typename outlier_finder_t>
struct Chi2FitterOptions {
  using Calibrator = calibrator_t;
  using OutlierFinder = outlier_finder_t;

  /// PropagatorOptions with context.
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param calibrator_t The source link calibrator
  /// @param outlierFinder_ The outlier finder
  /// @param logger_ The logger wrapper
  /// @param pOPtions The plain propagator options
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  Chi2FitterOptions(const GeometryContext& gctx,
                      const MagneticFieldContext& mctx,
                      std::reference_wrapper<const CalibrationContext> cctx,
                      Calibrator calibrator_, OutlierFinder outlierFinder_,
                      LoggerWrapper logger_,
                      const PropagatorPlainOptions& pOptions,
                      const Surface* rSurface = nullptr,
                      bool mScattering = false, bool eLoss = false)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        calibrator(std::move(calibrator_)),
        outlierFinder(std::move(outlierFinder_)),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        logger(logger_) {}
  /// Contexts are required and the options must not be default-constructible.
  Chi2FitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The source link calibrator.
  Calibrator calibrator;

  /// The outlier finder.
  OutlierFinder outlierFinder;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = true;

  /// Whether to consider energy loss
  bool energyLoss = true;

  /// Logger
  LoggerWrapper logger;
};

template <typename source_link_t>
struct Chi2FitterResult {

  // counter for the handled states (incremented inside processSurface())
  size_t processedStates = 0;

  // Indicator if track fitting has been done
  bool finished = false;

  // count the number of traversed planes = $m$
  size_t planeCounter = 0;




  // r = M - h(X)

  // X is the vector of starting track parameters, 6x1. Also has an associated 6x6 cov matrix.
  ActsMatrix<6,1> X;

  // H is the measurement function, which projects the track parameters to the
  // measurement. 2nx6 matrix.
  ActsDynamicMatrix H;

  // M is a 2xn matrix,  (2nx1??)
  ActsDynamicMatrix M;

  Result<void> result{Result<void>::success()};

  // ActsDynamicMatrix fullGlobalTrackParamsCov(nSmoothedStates* eBoundSize,
  //                                            nSmoothedStates* eBoundSize);
  // fullGlobalTrackParamsCov.setZero();
};


/// Chi square fitter implementation.
template <typename propagator_t>
class Chi2Fitter {
  /// The navigator type
  using Chi2Navigator = typename propagator_t::Navigator;

 public:
  Chi2Fitter(propagator_t pPropagator)
      : m_propagator(std::move(pPropagator)) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// @brief Propagator Actor plugin for the Chi2Filter
  ///
  /// @tparam source_link_t is an type fulfilling the @c SourceLinkConcept
  /// @tparam parameters_t The type of parameters used for "local" paremeters.
  /// @tparam calibrator_t The type of calibrator
  /// @tparam outlier_finder_t Type of the outlier finder class
  ///
  /// The Chi2Actor does not rely on the measurements to be
  /// sorted along the track.
  template <typename source_link_t, typename parameters_t,
            typename calibrator_t, typename outlier_finder_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = Chi2FitterResult<source_link_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    std::map<GeometryIdentifier, source_link_t> inputMeasurements;

    /// Whether to consider multiple scattering.
    bool multipleScattering = false; // TODO: add later

    /// Whether to consider energy loss.
    bool energyLoss = false; // TODO: add later

    /// The Surface beeing
    SurfaceReached targetReached;



    /// @brief Chi square actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
      const auto& logger = state.options.logger;

      if (result.finished) {
        return;
      }

      ACTS_VERBOSE("Chi2Fitter step");

      result.planeCounter += 1;

      // Update:
      // - Waiting for a current surface
      auto surface = state.navigation.currentSurface;
      if (surface != nullptr) {
          ACTS_VERBOSE("SURFACE: process chi2...");
          auto res = processSurface(surface, state, stepper, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in  filter: " << res.error());
             result.result = res.error();
           }
      } else {
        ACTS_VERBOSE("not on surface...")
      }

      // TODO: set `result.finished = true` after verifying if targetReached(state, stepper, *targetSurface)


    }

    /// @brief Chi2 actor operation : process surface
    ///
    /// (corresponds to KalmanFitter's `filter` method.)
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> processSurface(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, result_type& result) const {
      const auto& logger = state.options.logger;
      
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = inputMeasurements.find(surface->geometryId());
      if (sourcelink_it != inputMeasurements.end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");

        // Transport the covariance to the surface
        stepper.covarianceTransport(state.stepping, *surface);

        // Bind the transported state to the current surface
        auto [boundParams, jacobian, pathLength] =
            stepper.boundState(state.stepping, *surface, false);


        // -->>> here, the KalmanUpdate is run. <<<--


        // We count the processed state
        ++result.processedStates;
      }
      
      return Result<void>::success();
    }

  };

  template <typename source_link_t, typename parameters_t,
            typename calibrator_t, typename outlier_finder_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type =
        Actor<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      if (!result.result.ok() or result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation of the foward filter, calls the
  /// the filter and smoother/reversed filter
  ///
  /// @tparam source_link_t Type of the source link
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam calibrator_t Type of the source link calibrator
  /// @tparam outlier_finder_t Type of the outlier finder
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param chi2FitterOptions Chi2FitterOptions steering the fit
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_t, typename start_parameters_t,
            typename calibrator_t, typename outlier_finder_t,
            typename parameters_t = BoundTrackParameters>
  auto fit(const std::vector<source_link_t>& sourcelinks,
           const start_parameters_t& sParameters,
           const Chi2FitterOptions<calibrator_t, outlier_finder_t>& chi2FitterOptions) {
    const auto& logger = chi2FitterOptions.logger;

    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::map<GeometryIdentifier, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      inputMeasurements.emplace(sl.geometryId(), sl);
    }





    // Create the ActionList and AbortList
    using Chi2Aborter =
        Aborter<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;
    using Chi2Actor =
        Actor<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;
    using Chi2Result = typename Chi2Actor::result_type;
    using Actors = ActionList<Chi2Actor>;
    using Aborters = AbortList<Chi2Aborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> propOptions(
        chi2FitterOptions.geoContext, chi2FitterOptions.magFieldContext, logger);

    // Set the trivial propagator options
    propOptions.setPlainOptions(chi2FitterOptions.propagatorPlainOptions);

    // Catch the actor and set the measurements
    auto& chi2Actor = propOptions.actionList.template get<Chi2Actor>();
    chi2Actor.inputMeasurements = std::move(inputMeasurements);
    chi2Actor.targetSurface = chi2FitterOptions.referenceSurface;
    chi2Actor.multipleScattering = chi2FitterOptions.multipleScattering;
    chi2Actor.energyLoss = chi2FitterOptions.energyLoss;

    // chi2Actor.m_calibrator = chi2FitterOptions.calibrator;
    // chi2Actor.m_outlierFinder = chi2FitterOptions.outlierFinder;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, propOptions);

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error());
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto chi2Result = propRes.template get<Chi2Result>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (chi2Result.result.ok() and not chi2Result.processedStates) {
      chi2Result.result = Result<void>(KalmanFitterError::NoMeasurementFound); // TODO: generic Fitter error?
    }

    if (!chi2Result.result.ok()) {
      ACTS_ERROR("Chi2Filter failed: " << chi2Result.result.error());
      return chi2Result.result.error();
    }

    // Return the converted Track
    return chi2Result;
  }


};

}  // namespace Acts
