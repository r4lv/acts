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
  bool multipleScattering = false;  // TODO: add later

  /// Whether to consider energy loss
  bool energyLoss = false;  // TODO: add later

  /// Logger
  LoggerWrapper logger;
};

template <typename source_link_t>
struct Chi2FitterResult {

  MultiTrajectory<source_link_t> fittedStates; // workaround. chi2 should use another trajectory type

  // counter for the handled states (incremented inside processSurface())
  size_t processedStates = 0;

  // Indicator if track fitting has been done
  bool finished = false;



  // r = M - h(X)

  // (1)
  ActsDynamicVector r_measurement_residual; // nx1
  ActsDynamicMatrix rcov_measurement_covariance; // sigma_meas

  std::vector <Acts::ActsVector<2>> m_measurements_per_surface; // here: 2n x 1 (or nx2 better?)

          // (2)
          // Matrix: derivative of residuals wrt starting parameters
          // Vector: ∂ χ²  wrt starting track parameters
          // Matrix: ∂² χ² wrt starting track parameters

          // (3)
          // Vector: update of starting track parameters

          // ————————————————————————————————————————————————————————————

          // X is the vector of starting track parameters, 6x1. Also has an
          // associated 6x6 cov matrix.
  // ActsVector<6> X_track_parameters;
  // ActsMatrix<6,6> Xcov_track_parameters_cov;

  // H is the measurement function, which projects the track parameters to the
  // measurement. 2nx6 matrix.
  // ActsDynamicMatrix H_measurement_function;

  // M is a 2nx1 matrix,  (why 2nx1 and not nx2?)
  // m measurements, each is 2dim, with a 2x2 cov matrix
  ActsDynamicVector M_measurement;
  ActsDynamicMatrix Mcov_measurement_cov; // 2mx2?
  // std::vector<Vector3> position;  ?
  

  // ActsDynamicMatrix fullGlobalTrackParamsCov(nSmoothedStates* eBoundSize,
  // nSmoothedStates* eBoundSize); fullGlobalTrackParamsCov.setZero();

  // ————————————————————————————————————————————————————————————

  Result<void> result{Result<void>::success()};


};


/// Chi2 fitter implementation.
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

      // Add the measurement surface as external surface to navigator.
      // We will try to hit those surface by ignoring boundary checks.
      if (result.processedStates == 0) {
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
          auto res = processSurface(surface, state, stepper, result);
          if (!res.ok()) {
            ACTS_ERROR("      Error in processSurface: " << res.error());
             result.result = res.error();
           }
      } else {
        ACTS_INFO("   not on surface. Done")
      }

      if (targetSurface == nullptr){
        // in KF: "no target surface provided: fitting is finished here" WHY?
        // result.finished = true;
      }
      else if (targetReached(state, stepper, *targetSurface)){
        ACTS_INFO("operator() reached target surface. Completed fitting.");

        // // Transport & bind the parameter to the final surface
        // auto fittedState = stepper.boundState(state.stepping, *targetSurface);
        // // Assign the fitted parameters
        // result.fittedParameters = std::get<BoundTrackParameters>(fittedState);
      }

      // ACTS_INFO("finishing Actor operator(): result.processedStates = " << result.processedStates);

      if (result.processedStates == 6){
        // TODO: this is a quick hack to stop the propagator.
        // TODO: use targetReached(state, stepper, *targetSurface) like in KF
        ACTS_INFO("operator(): 6 surfaces processed -> finishing");
        result.finished = true;
      }
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
      
      // Try to find the surface in all measurement surfaces
      auto sourcelink_it = inputMeasurements.find(surface->geometryId());
      // inputMeasurements is a std::map<GeometryIdentifier, source_link_t>
      if (sourcelink_it != inputMeasurements.end()) {
        ++result.processedStates;
        
        ACTS_INFO("processSurface() Measurement surface " << surface->geometryId() << " detected.");

        // Transport the covariance to the surface
        stepper.covarianceTransport(state.stepping, *surface);

        // Bind the transported state to the current surface
        auto [boundParams, jacobian, pathLength] = stepper.boundState(state.stepping, *surface, false);
        // EigenStepper.boundState() returns a BoundState,
        // which is a tuple<Acts::BoundTrackParameters, Acts::BoundMatrix, double>
        // Acts:BoundTrackParameters = Acts::SingleBoundTrackParameters<Acts::SinglyCharged>
        //            from Acts/EventData/SingleBoundTrackParameters.hpp
        // Acts::BoundMatrix = Acts::ActsMatrix<6U, 6U>


        // TODO: use `auto` everywhere. Concrete types are kept for debugging

        const Acts::BoundVector& predictedParams = std::move(boundParams.parameters());
        // `Acts::BoundVector` is a `Acts::ActsVector<6U>`
        const Acts::BoundSymMatrix predictedCovariance = std::move(*boundParams.covariance()); // TODO: why move?

        // Jacobian between measurement A and measurement B

        ACTS_INFO("   collecting information");
        ACTS_INFO("      parameters: " << predictedParams.transpose());
        // ACTS_INFO("      jacobian: " << jacobian);
        // ACTS_INFO("      pathLength: " << pathLength);

        Acts::Test::TestSourceLink uncalibrated = sourcelink_it->second;
        // source_link_t is Acts::Test::TestSourceLink

        Acts::ActsVector<2> measurementParameters = uncalibrated.parameters;
        Acts::ActsSymMatrix<2> measurementCovariance = uncalibrated.covariance;
        constexpr size_t kMeasurementSize = decltype(measurementParameters)::RowsAtCompileTime;

        // TestSourceLinkCalibrator
        // the 

        // TODO: how to store the calibrated measurements? custom trackStateProxy?
        //
        // from KF:
        // We have predicted parameters, so calibrate the uncalibrated input measuerement
        // Acts::Measurement calibrated;
        ACTS_INFO("      calibration visit...");



        
        

        const auto& calibratorResult = m_calibrator(uncalibrated, predictedParams);
        // returns a BoundVariantMeasurement<TestSourceLink>

        // using Meas_t = Acts::Measurement<Acts::Test::TestSourceLink, Acts::BoundIndices, kMeasurementSize>;
        // const auto& calibrated = std::get<Acts::Measurement<Acts::Test::TestSourceLink, Acts::BoundIndices, 6>>(calibratorResult);
        // const auto& calibrated = std::get<1>(calibratorResult);
        // const auto& calibrated = std::get<Acts::Measurement>(calibratorResult);
        // const Acts::Measurement<Acts::Test::TestSourceLink, Acts::BoundIndices, 6> calibrated

        // static constexpr size_t kSize = 6;
        // static constexpr size_t kFullSize = 6;
        // using ProjectionMatrix = ActsMatrix<kSize, kFullSize>;

        using Covariance =  Acts::detail_lt::Types<eBoundSize, true>::CovarianceMap;
        // without typename!

        const size_t M = 6; // maximum number of measurement dimensions. TODO: where can I get this from?
        constexpr static auto ProjectorFlags = Eigen::RowMajor | Eigen::AutoAlign;
        using Projector = Eigen::Matrix<typename Covariance::Scalar, M, eBoundSize, ProjectorFlags>;


        // /scratch/rfarkas/phd/Code/acts-chisquare/spack-env-4/.spack-env/view/include/eigen3
        // /Eigen/src/Core/AssignEvaluator.h: In instantiation of 'void Eigen::internal::call_
        // assignment_no_alias(Dst&, const Src&, const Func&) [with Dst = Eigen::Matrix<double
        // , 6, 6, 1>; Src = Eigen::Matrix<double, 1, 6, 1, 1, 6>; Func = Eigen::internal::ass
        // ign_op<double, double>]':
        // /scratch/rfarkas/phd/Code/acts-chisquare/spack-env-4/.spack-env/view/include/eigen3
        // /Eigen/src/Core/PlainObjectBase.h:732:41:   required from 'Derived& Eigen::PlainObj
        // ectBase<Derived>::_set_noalias(const Eigen::DenseBase<OtherDerived>&) [with OtherDe
        // rived = Eigen::Matrix<double, 1, 6, 1, 1, 6>; Derived = Eigen::Matrix<double, 6, 6,
        // 1>]'

        std::visit(
            [&](const auto& calibrated_meas) {
              // 
              // calibrated_meas is of type Acts::Measurement<Acts::Test::TestSourceLink, Acts::BoundIndices, 6>
              const auto& proj = calibrated_meas.projector();
              // const Projector& proj2 = calibrated_meas.projector();

              Projector fullProjector = decltype(fullProjector)::Zero();
              // fullProjector.template topLeftCorner<rows, cols>() = proj;
              // error: 'RowsAtCompileTime' is not a member of 'const Eigen::Matrix<double, 1, 6, 1, 1, 6>&'
              
              // Acts::ActsVector<2> measurementParameters = uncalibrated.parameters;
              // Acts::ActsSymMatrix<2> measurementCovariance = uncalibrated.covariance;
              // constexpr size_t kMeasurementSize = decltype(measurementParameters)::RowsAtCompileTime;

              constexpr int rows = decltype(proj)::RowsAtCompileTime;
              // RowsAtCompileTime' is not a member of 'const Eigen::Matrix<double, 1, 6, 1, 1, 6>&'
              constexpr int cols = decltype(proj)::ColsAtCompileTime;

              // proj is `const Eigen::MatrixBase<Derived>&`
              ACTS_INFO("         got calibrated measurement + projector");
              ACTS_INFO("         rows=" << rows << "  cols=" << cols);
              // ACTS_INFO("         projector " << decltype(proj));
              ACTS_INFO("         projector " << proj);
              ACTS_INFO("         projector 2 " << fullProjector);
              // ACTS_INFO("         matrix rows: " << decltype(proj)::RowsAtCompileTime);
              // ACTS_INFO("         matrix cols: " << decltype(proj)::ColsAtCompileTime);
              // trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(uncalibrated, predictedParams));
        ACTS_INFO("      done calibrating.");

        


        ACTS_INFO("      measurement parameters:" << measurementParameters.transpose());
        result.m_measurements_per_surface.push_back(measurementParameters);
        ACTS_INFO("      measurement dimension:" << kMeasurementSize);

        // the Acts::Measurement<SourceLink, BoundIndices, kMeasurementSize>& meas object has a .projector() method

        // from GainMatrixUpdater.hpp:78
        //
        // const auto H =
        //     trackState.projector()
        //         .template topLeftCorner<kMeasurementSize, eBoundSize>()
        //         .eval();
        // 
        // ACTS_INFO("      Measurement projector H:\n" << H);
        // 
        // const auto K =
        //     (predictedCovariance * H.transpose() *
        //      (H * predictedCovariance * H.transpose() + calibratedCovariance)
        //          .inverse())
        //         .eval();
        
        // GainMatrixUpdater.hpp :118
        //
        // ParametersVector residual = calibrated - H * filtered;
        // [...]
        // trackState.chi2() =
        //     (residual.transpose() *
        //      ((CovarianceMatrix::Identity() - H * K) * calibratedCovariance)
        //          .inverse() *
        //      residual)
        //         .value();




// 
        // auto residuals = state.calibrated() - state.projector() *  state.predicted();

        // TODO: collect measurement information
        //   - residual vector
        //  (r^2_meas, sigma^2_meas)
        
        // TODO: calculate residual + chisq, and derivatives of this measurement
        // TODO: calculate update of track parameters



        // -->>> in the Kalman Fitter, we run the KalmanUpdate m_updater() here <<<--
      }
      
      ACTS_INFO("processSurface(): returning Result::success()");
      return Result<void>::success();
    }

    /// The measurement calibrator
    calibrator_t m_calibrator;
    // here, calibrator_t is a TestSourceLinkCalibrator

    /// The outlier finder
    outlier_finder_t m_outlierFinder;
  };

  template <typename source_link_t, typename parameters_t,
            typename calibrator_t, typename outlier_finder_t>
  class Aborter {
   public:
    /// Broadcast the action_type
    using action_type =
        Actor<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      // TODO: we have no logger in the aborter
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
           const Chi2FitterOptions<calibrator_t, outlier_finder_t>& chi2FitterOptions) const
      -> Result<Chi2FitterResult<source_link_t>> {
    const auto& logger = chi2FitterOptions.logger;

    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_INFO("Preparing " << sourcelinks.size() << " input measurements");
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
    chi2Actor.multipleScattering = chi2FitterOptions.multipleScattering;
    chi2Actor.energyLoss = chi2FitterOptions.energyLoss;

    chi2Actor.m_calibrator = chi2FitterOptions.calibrator;
    chi2Actor.m_outlierFinder = chi2FitterOptions.outlierFinder;

    // start propagator with initial estimte, let actor collect necessary parameters
    auto result = m_propagator.template propagate(sParameters, propOptions);
    ACTS_INFO("propagation done");
    
    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error()); // PropagatorError:3
      return result.error();
    }


    // fitter retreives the object
    // runs minimizer (ideally we can make it pluggable)
    // if it converges, fine. Else, the fitter needs to re-invoke the proapgator from the new start and go ahead.


    /// Get the result of the fit
    const auto& propRes = *result;
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

    ACTS_INFO("results.m_measurements_per_surface contains " << chi2Result.m_measurements_per_surface.size() << " results.");
    if (chi2Result.m_measurements_per_surface.size() > 1){
      ACTS_INFO("   surf0: " << chi2Result.m_measurements_per_surface[0].transpose());
      ACTS_INFO("   surf1: " << chi2Result.m_measurements_per_surface[1].transpose());
    }
    

    // Return the converted track
    return chi2Result;

  }
};

}  // namespace Acts
