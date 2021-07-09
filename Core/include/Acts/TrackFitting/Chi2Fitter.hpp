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
// TODO: do we need a Chi2FitterError, or would a more generic FitterError be sufficient?
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
  /// @param calibrator_ The source link calibrator
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
                    const Surface* rSurface = nullptr, bool mScattering = false,
                    bool eLoss = false, int nIter = 1,
                    bool calcFinalChi2_ = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        calibrator(std::move(calibrator_)),
        outlierFinder(std::move(outlierFinder_)),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        nUpdates(nIter),
        calcFinalChi2(calcFinalChi2_),
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
  bool multipleScattering = false;  // TODO: add later

  /// Whether to consider energy loss
  bool energyLoss = false;  // TODO: add later

  /// Number of iterations to improve chi2
  int nUpdates;

  /// Whether to do a last propagation step, just to get the latest chi2 value
  bool calcFinalChi2;

  /// Logger
  LoggerWrapper logger;
};

template <typename source_link_t>
struct Chi2FitterResult {

  MultiTrajectory<source_link_t> fittedStates;

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t trackTip = SIZE_MAX;

  // The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  // Counter for states with measurements
  size_t measurementStates = 0; // TODO: use to handle finalization

  // counter for the handled states (incremented inside processSurface())
  size_t processedStates = 0;


  // Indicator if track fitting has been done
  bool finished = false;
  // we do not have smoothed, reset, reversed, like the KF does.

  // Measurement surfaces without hits
  // std::vector<const Surface*> missedActiveSurfaces; // not used for now

  // Measurement surfaces handled in both forward and backward filtering
  // std::vector<const Surface*> passedAgainSurfaces; // not used for now

  // intermediary results / collectors for chi2 calculation:

  // std::vector<ActsVector<2>> r_measurement_residuals;
  // std::vector<ActsVector<2>> rcov_measurement_covariance;
  // std::vector<BoundVector> H_projector_collector;
  // ActsDynamicVector M_measurement;
  // ActsDynamicMatrix Mcov_measurement_cov; // 2mx2?

  // collectors

  std::vector<ActsScalar> m_measurements_collector; // pixel layers add 2 entries, strip layers only one.
  std::vector<ActsScalar> mcov_measurement_covariance_collector; // we assume measurements are not correlated
  std::vector<ActsScalar> r_residual_collector;

  BoundVector d1_deriv_sum_collector = BoundVector::Zero(); // first derivative of chi2 wrt starting track parameters, 6x1
  BoundMatrix d2_deriv_sum_collector = BoundMatrix::Zero();

  BoundMatrix jacobianFromStart = BoundMatrix::Identity();

  // chi2 results
  ActsDynamicVector residuals;
  ActsDynamicMatrix covariance;
  ActsScalar chisquare = -1;

  std::vector<ActsScalar> chisquares;

  Result<void> result{Result<void>::success()};
};


/// Chi2 fitter implementation.
template <typename propagator_t>
class Chi2Fitter {
  using Chi2Navigator = typename propagator_t::Navigator;

 public:
  Chi2Fitter(propagator_t pPropagator)
      : m_propagator(std::move(pPropagator)) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// @brief Propagator Actor plugin for the Chi2Filter
  ///
  /// @tparam source_link_t is a type fulfilling the @c SourceLinkConcept
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

    /// The measurement calibrator
    calibrator_t m_calibrator;

    /// The outlier finder
    outlier_finder_t m_outlierFinder;



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
        ACTS_INFO("CHI2 operator() result.finished...");
        return;
      }

      ACTS_VERBOSE("Chi2Fitter step");

      // Add the measurement surface as external surface to navigator.
      // We will try to hit those surface by ignoring boundary checks.
      if (result.processedStates == 0) {
        ACTS_INFO("CHI2 operator() no processed states yet...");
        for (auto measurementIt = inputMeasurements.begin();
             measurementIt != inputMeasurements.end(); measurementIt++) {
          state.navigation.externalSurfaces.insert(
              std::pair<uint64_t, GeometryIdentifier>(
                  measurementIt->first.layer(), measurementIt->first));
        }
      }

      // wait for surface
      auto surface = state.navigation.currentSurface;
      if (surface != nullptr) {
        ACTS_INFO("CHI2 operator() we have a surface, calling processSurface with surface " << surface);
        auto res = processSurface(surface, state, stepper, result);
        if (!res.ok()) {
          ACTS_ERROR("      Error in processSurface: " << res.error());
          result.result = res.error();
           }
      } else {
        ACTS_INFO("   not on surface. Done");
      }

      ACTS_INFO("CHI2 operator() after processSurface()");

      if (targetSurface == nullptr){
        // in KF: "no target surface provided: fitting is finished here" WHY?
        // result.finished = true;
      }
      // else if (targetReached(state, stepper, *targetSurface)){
      //   ACTS_INFO("operator() reached target surface. Completed fitting. TODO");
      // 
      //   // // Transport & bind the parameter to the final surface
      //   // auto fittedState = stepper.boundState(state.stepping, *targetSurface);
      //   // // Assign the fitted parameters
      //   // result.fittedParameters = std::get<BoundTrackParameters>(fittedState);
      // }

      ACTS_INFO("CHI2 after targetReached");

      // ACTS_INFO("finishing Actor operator(): result.processedStates = " << result.processedStates);

      // if (result.processedStates == 6){
      //   // TODO: this is a quick hack to stop the propagator for the Unit Test.
      //   // TODO: use targetReached(state, stepper, *targetSurface) like in KF
      //   ACTS_INFO("Actor.operator(): 6 surfaces processed -> finishing FIXME");
      //   result.finished = true;
      // }
    }

    /// @brief Chi2 actor operation: process surface
    ///
    /// (corresponds to KalmanFitter's `filter` method.)
    ///
    /// @tparam propagator_state_t is the type of Propagator state
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

      ACTS_VERBOSE("processSurface start");
      //  << surface->geometryId());

      // Try to find the surface in all measurement surfaces
      auto sourcelink_it = inputMeasurements.find(surface->geometryId());
      // inputMeasurements is a std::map<GeometryIdentifier, source_link_t>
      if (sourcelink_it != inputMeasurements.end()) {
        ++result.processedStates;

        ACTS_VERBOSE("  Measurement surface " << surface->geometryId()
                                            << " detected.");

        // Transport the covariance to the surface
        stepper.covarianceTransport(state.stepping, *surface);

        // Update state and stepper with pre material effects
        // materialInteractor(surface, state, stepper, preUpdate);

        // Bind the transported state to the current surface
        auto [boundParams, jacobian, pathLength] = stepper.boundState(state.stepping, *surface, false); // last parameter: do not transport cov
        // EigenStepper.boundState() returns a tuple<BoundTrackParameters, BoundMatrix, double>
        // Acts:BoundTrackParameters = Acts::SingleBoundTrackParameters<Acts::SinglyCharged>
        //            from Acts/EventData/SingleBoundTrackParameters.hpp
        // BoundMatrix is Acts::ActsMatrix<6U, 6U> (Jacobian between measurement A and measurement B)
        // here, the jacobian is "the stepwise jacobian towards the bound state (from last bound)".
        // We need the full jacobian, from the first surface to the current.

        



        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)
        result.trackTip = result.fittedStates.addTrackState(
            TrackStatePropMask::All, result.trackTip);

        // now get track state proxy back
        auto trackStateProxy =
            result.fittedStates.getTrackState(result.trackTip);

        trackStateProxy.setReferenceSurface(surface->getSharedPtr());

        // assign the source link to the track state
        trackStateProxy.uncalibrated() = sourcelink_it->second;

        // Fill the track state
        trackStateProxy.predicted() = std::move(boundParams.parameters());
        trackStateProxy.predictedCovariance() =
            std::move(*boundParams.covariance());
        trackStateProxy.jacobian() = std::move(jacobian);
        trackStateProxy.pathLength() = std::move(pathLength);



        result.jacobianFromStart = jacobian * result.jacobianFromStart;
        // TODO: also update result.jacobianFromStart when we have no measurement at current surface?
        

        // const BoundVector& predictedParams = std::move(boundParams.parameters());

        // auto uncalibrated = sourcelink_it->second;
        // ^ in the Unit Test, this is a Acts::Test::TestSourceLink

        // ActsVector<2> measurementParameters = uncalibrated.parameters;
        // Acts::ActsSymMatrix<2> measurementCovariance = uncalibrated.covariance;

        // calibrate the uncalibrated input measurement + collect information
        std::visit(
            [&](const auto& calibrated) {

              const auto& proj = calibrated.projector(); // simple projection matrix H_is, composed of 1 and 0
              // proj is a Eigen::MatrixBase<Derived>. here, proj can be 1x6 (strip) or 2x6 (pixel)


              const auto& Hi = proj * result.jacobianFromStart; // has dimension 1x6 or 2x6

              const auto& parameters = calibrated.parameters();
              // ParametersVector = ActsVector

              const auto& covariance = calibrated.covariance(); // 2x2 or 1x1. Should be diagonal.
              const auto& covarianceInverse = covariance.inverse();
              // const CovarianceMatrix& = ActSymMatrix

              const auto& residuals = calibrated.residuals(
                  trackStateProxy.predicted());  // 2x1 or 1x1
              //  ParametersVector residuals(const FullParametersVector&
              //  reference)

              const auto& deriv1 = -2 * Hi.transpose() * covarianceInverse * residuals; // 6x1
              const auto& deriv2 = 2*Hi.transpose() * covarianceInverse * Hi; // 6x6

              result.d1_deriv_sum_collector += deriv1;
              result.d2_deriv_sum_collector += deriv2;

              result.m_measurements_collector.push_back(parameters(0));
              result.mcov_measurement_covariance_collector.push_back(covariance(0,0));

              result.r_residual_collector.push_back(residuals(0));
              // predictedParams
              // ACTS_INFO("         calculating residuals. Predicted Params:\n" <<predictedParams << "\n----\n" << parameters);

              // ACTS_INFO("  residuals: " << residuals.transpose());

              trackStateProxy.setCalibrated(calibrated);

              if (parameters.rows() == 2){
                result.m_measurements_collector.push_back(parameters(1));
                // This however depends on the matrix's storage order. All Eigen matrices default to column-major storage order, but this can be changed to row-major, see Storage orders.
                result.r_residual_collector.push_back(residuals(1));
                result.mcov_measurement_covariance_collector.push_back(covariance(0, 0));
              }
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.predicted()));

 
        // Get and set the type flags
        auto& typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::MaterialFlag);
        typeFlags.set(TrackStateFlag::ParameterFlag);


        if (true){
          // TODO: MeasurementFlag should only be check if we don't have an
          // outlier. the KF checks it with `m_outlierFinder(trackStateProxy)`.
          typeFlags.set(TrackStateFlag::MeasurementFlag);
        } else {
          ACTS_VERBOSE("Measurement is determined to be an outlier.");
          typeFlags.set(TrackStateFlag::OutlierFlag);
        }

        


      } else if (surface->surfaceMaterial() != nullptr) {
        ACTS_INFO("TODO: hole / passive material"); // TODO: see KF
      } else {
        ACTS_INFO("CHI2 processSurface else! TODO not handled!");
      }

      return Result<void>::success();
    }

  };

  template <typename source_link_t, typename parameters_t,
            typename calibrator_t, typename outlier_finder_t>
  // DEBUG: for the unit test, we use Acts::Test::TestSourceLink as source_link_t
  class Aborter {
   public:
    /// Broadcast the action_type
    using action_type =
        Actor<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& state, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      const auto& logger = state.options.logger;
      if (!result.result.ok() or result.finished) {
        ACTS_INFO("CHI2 Aborter: return true");
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
           const Chi2FitterOptions<calibrator_t, outlier_finder_t>& chi2FitterOptions) const
      -> Result<Chi2FitterResult<source_link_t>> {
    const auto& logger = chi2FitterOptions.logger;

    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_INFO("CHI2 Preparing " << sourcelinks.size() << " input measurements");
    std::map<GeometryIdentifier, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      inputMeasurements.emplace(sl.geometryId(), sl);
    }

    // size_t n = sourcelinks.size(); // required for settings the size of the dynamic result matrices
    // TODO: for now, we use STL objects to collect the information during propagation.

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
    // ^ what is this?
    chi2Actor.multipleScattering = chi2FitterOptions.multipleScattering;
    chi2Actor.energyLoss = chi2FitterOptions.energyLoss;
    // KF: reversedFiltering = kfOptions.reversedFiltering;
    chi2Actor.m_calibrator = chi2FitterOptions.calibrator;
    chi2Actor.m_outlierFinder = chi2FitterOptions.outlierFinder;

    // overwrite
    // chi2FitterOptions.calcFinalChi2 = true;
    // chi2FitterOptions.nUpdates = 2;

    // std::variant<start_parameters_t, parameters_t> vParams = sParameters;
    using paramType = typename std::conditional<
        std::is_same<start_parameters_t, parameters_t>::value,
        std::variant<start_parameters_t>,
        std::variant<start_parameters_t, parameters_t>>::type;
    // std::variant<start_parameters_t> vParams = sParameters;
    paramType vParams = sParameters;
    // TODO: in the example, start_parameters_t and parameters_t are the same!!!!


    Chi2Result c2r; // the result object which will be returned


    // nUpdates = 1 (default) -> [update, chi2, delta] + [update, chi2] -> 3 propagations
    // nUpdates = 2 -> [up, chi2, ∆] + [up, chi2, ∆] + [update, chi2] -> 3 propagations

    ACTS_INFO("CHI2 running PRE propagator...");
    auto resultPre = m_propagator.template propagate(sParameters, propOptions);
    ACTS_INFO("CHI2 PRE propagator done.");



    for (int i = 0; i <= chi2FitterOptions.nUpdates; ++i) {
      ACTS_INFO("CHI2 running propagator for iteration i=" << i << "...");
     
      auto result = std::visit(
          [mprop = m_propagator, propOptions](auto&& arg) {
            auto result_ = mprop.template propagate(arg, propOptions);
            // TODO: verify result.ok(), and return Chi2 result directly?
            return result_;
          },
          vParams);
      ACTS_VERBOSE("CHI2 finished propagation for iteration i=" << i <<".");

      if (!result.ok()) {
        ACTS_ERROR("Propapation failed: " << result.error());
        return result.error();
      }
      c2r = result.value().template get<Chi2Result>();

      /// It could happen that the fit ends in zero processed states.
      /// The result gets meaningless so such case is regarded as fit failure.
      if (c2r.result.ok() and not c2r.processedStates) {
        c2r.result = Result<void>(KalmanFitterError::NoMeasurementFound);
        // TODO: generic Fitter error instead of using the KalmanFitterError?
      }

      if (!c2r.result.ok()) {
        ACTS_ERROR("Chi2Filter failed: " << c2r.result.error());
        return c2r.result.error();
      }


      c2r.residuals =
          Eigen::Map<ActsDynamicVector>(c2r.r_residual_collector.data(),
                                        c2r.r_residual_collector.size());
      // TODO: is this safe? from Stackoverflow: "dangerous! Because the Eigen object will NOT create its own memory. It will operate on the memory provided by
      // "data". In other words, working with the Eigen object when the "data" object is out of scope will result in a segmentation fault (or memory access violation)."

      ActsDynamicVector variance = Eigen::Map<ActsDynamicVector>(
          c2r.mcov_measurement_covariance_collector.data(),
          c2r.mcov_measurement_covariance_collector.size());
      c2r.covariance = variance.asDiagonal();
      // this should match the `MeasurementResolution`s from the test

      // // build Matrix from projector list
      // ActsDynamicMatrix H = ActsDynamicMatrix::Zero(n, 6); // TODO: as we overwrite the values, use another constructor?
      // // error: static assertion failed: YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE _TO_CAST_NUMERIC_TYPES_EXPLICITLY
      // for(int i = 0; i < n; ++i) {
      //   ACTS_INFO("       H projector i="<<i<<" : " << c2r.H_projector_collector[i].transpose());
      //   H.row(i) = c2r.H_projector_collector[i];
      // }
      

      c2r.chisquare = c2r.residuals.transpose() *
                            c2r.covariance.inverse() *
                            c2r.residuals;
      c2r.chisquares.push_back(c2r.chisquare);
      ACTS_INFO("χ² = " << c2r.chisquare);

      if (i==chi2FitterOptions.nUpdates)
        break; // don't update parameters in last iteration

      // calculate delta on parameters
      BoundVector delta_start_parameters =
          c2r.d2_deriv_sum_collector.colPivHouseholderQr().solve(
              c2r.d1_deriv_sum_collector);

      c2r.fittedParameters = std::visit(
          [delta_start_parameters, logger](auto&& prevParams) {
            BoundVector newParamsVec =
                prevParams.parameters() - delta_start_parameters;
            ACTS_INFO("∆parameters = " << delta_start_parameters.transpose());
            ACTS_INFO("  updated parameters: " << newParamsVec.transpose());

            return BoundTrackParameters(
                prevParams.referenceSurface().getSharedPtr(),
                newParamsVec, prevParams.covariance());
          },
          vParams);

      vParams = c2r.fittedParameters.value(); // passed to next iteration

      // TODO: we need to add the fitted parameters to the TrackStateProxy, as 'smoothed'
      // so that the performance writer can use it.

    }  // <- end for loop

    // TODO: what should we do with the MultiTrajectory? Where do I set the 'smoothed'?

    ACTS_INFO("finished chi2 calculation, returning result...")

    // Return the converted track
    return c2r;
  }


};

}  // namespace Acts
