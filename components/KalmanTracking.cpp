#include <iostream>
#include "KalmanTracking.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/SurfaceHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/MaterialManager.h"
#include "podio/Frame.h"
#include "k4FWCore/PodioDataSvc.h"
#include "k4FWCore/DataHandle.h"

// Some physical constants
constexpr double electron_mass = 0.000511;  // GeV
constexpr double muon_mass = 0.10566;       // GeV
constexpr double pion_mass = 0.13957;       // GeV
constexpr double kaon_mass = 0.49368;       // GeV
constexpr double proton_mass = 0.93827;     // GeV
constexpr double c_light = 299792458.0;     // m/s
constexpr double qe = 1.602176634e-19;      // coulombs

//------------------------------------------------------------------------------
// TrackState Implementation
//------------------------------------------------------------------------------

TrackState::TrackState(const Eigen::Vector5d& params, const Eigen::Matrix5d& cov, 
                       const dd4hep::rec::Surface* surface)
    : _params(params), _cov(cov), _surface(surface) {
}

Eigen::Vector3d TrackState::momentum() const {
    // Extract track parameters
    double qOverPt = _params(0);  // q/pT
    double phi = _params(1);      // azimuthal angle
    double eta = _params(2);      // pseudorapidity
    
    // Convert eta to theta
    double theta = 2.0 * std::atan(std::exp(-eta));
    
    // Calculate pT, px, py, pz
    double pT = std::abs(1.0 / qOverPt);
    double px = pT * std::cos(phi);
    double py = pT * std::sin(phi);
    double pz = pT / std::tan(theta);
    
    return Eigen::Vector3d(px, py, pz);
}

Eigen::Vector3d TrackState::positionAtOrigin() const {
    // Extract track parameters
    double d0 = _params(3);      // transverse impact parameter
    double z0 = _params(4);      // longitudinal impact parameter
    double phi = _params(1);     // azimuthal angle
    
    // Calculate position at closest approach to origin
    double x0 = -d0 * std::sin(phi);
    double y0 = d0 * std::cos(phi);
    
    return Eigen::Vector3d(x0, y0, z0);
}

Eigen::Matrix<double, 2, 5> TrackState::projectionMatrix(const dd4hep::rec::Surface* surface) const {
    // Initialize projection matrix
    Eigen::Matrix<double, 2, 5> H = Eigen::Matrix<double, 2, 5>::Zero();
    
    // Get surface properties
    dd4hep::rec::Vector3D u_axis = surface->u();
    dd4hep::rec::Vector3D v_axis = surface->v();
    dd4hep::rec::Vector3D normal = surface->normal();
    dd4hep::rec::Vector3D origin = surface->origin(); // convert cm to mm
    
    // Extract track parameters
    double qOverPt = _params(0);
    double phi = _params(1);
    double eta = _params(2);
    double d0 = _params(3);
    double z0 = _params(4);
    
    // Calculate theta from eta
    double theta = 2.0 * std::atan(std::exp(-eta));
    
    // Calculate track direction
    Eigen::Vector3d dir;
    dir << std::cos(phi) * std::sin(theta),
           std::sin(phi) * std::sin(theta),
           std::cos(theta);
    
    // Calculate track reference position
    Eigen::Vector3d pos;
    pos << -d0 * std::sin(phi),
           d0 * std::cos(phi),
           z0;
    
    // Calculate the intersection point of the track with the surface
    // This is a simplified approximation
    dd4hep::rec::Vector3D normVec = normal.unit();
    dd4hep::rec::Vector3D surfOrigin = origin;
    double t = (normVec.dot(surfOrigin) - normVec.dot(dd4hep::rec::Vector3D(pos.x(), pos.y(), pos.z()))) / 
              normVec.dot(dd4hep::rec::Vector3D(dir.x(), dir.y(), dir.z()));
    
    // Check if intersection is valid
    if (t <= 0) {
        // If no forward intersection, use a simplified projection
        H(0, 1) = 10.0;  // Approximately project phi to u
        H(1, 2) = 10.0;  // Approximately project eta to v
        return H;
    }
    
    // Calculate intersection point
    Eigen::Vector3d intersectPos = pos + t * dir;
    
    // Calculate local coordinates of the intersection point
    dd4hep::rec::Vector3D globalPos(intersectPos.x(), intersectPos.y(), intersectPos.z());
    dd4hep::rec::Vector2D localPos = surface->globalToLocal(globalPos);
    
    // Now we need to calculate how changes in track parameters affect the local position
    // This requires computing partial derivatives
    
    // For a small change in each parameter, calculate the change in local position
    const double delta = 0.01;  // Small delta for numerical differentiation
    
    // For each track parameter, calculate partial derivatives
    for (int i = 0; i < 5; i++) {
        // Create a modified state with perturbed parameter
        Eigen::Vector5d modParams = _params;
        modParams(i) += delta;
        
        // Recalculate direction and position with perturbed parameter
        double modPhi = (i == 1) ? modParams(1) : phi;
        double modEta = (i == 2) ? modParams(2) : eta;
        double modD0 = (i == 3) ? modParams(3) : d0;
        double modZ0 = (i == 4) ? modParams(4) : z0;
        
        double modTheta = 2.0 * std::atan(std::exp(-modEta));
        
        Eigen::Vector3d modDir;
        modDir << std::cos(modPhi) * std::sin(modTheta),
                 std::sin(modPhi) * std::sin(modTheta),
                 std::cos(modTheta);
        
        Eigen::Vector3d modPos;
        modPos << -modD0 * std::sin(modPhi),
                 modD0 * std::cos(modPhi),
                 modZ0;
        
        // Calculate new intersection
        double modT = (normVec.dot(surfOrigin) - normVec.dot(dd4hep::rec::Vector3D(modPos.x(), modPos.y(), modPos.z()))) / 
                     normVec.dot(dd4hep::rec::Vector3D(modDir.x(), modDir.y(), modDir.z()));
        
        if (modT <= 0) continue;  // Skip if no valid intersection
        
        Eigen::Vector3d modIntersectPos = modPos + modT * modDir;
        dd4hep::rec::Vector3D modGlobalPos(modIntersectPos.x(), modIntersectPos.y(), modIntersectPos.z());
        dd4hep::rec::Vector2D modLocalPos = surface->globalToLocal(modGlobalPos);
        
        // Calculate derivatives (change in local position per change in parameter)
        H(0, i) = (modLocalPos.u() - localPos.u()) / delta;
        H(1, i) = (modLocalPos.v() - localPos.v()) / delta;
    }
    
    return H;
}

TrackState TrackState::predictTo(const dd4hep::rec::Surface* newSurface, 
                               const dd4hep::OverlayedField& field,
                               const dd4hep::rec::MaterialManager& matMgr,
                               const ParticleProperties& particle) const {
    try {
        // Extract current track parameters
        double qOverPt = _params(0);  // q/pT
        double phi = _params(1);      // azimuthal angle
        double eta = _params(2);      // pseudorapidity
        double d0 = _params(3);       // transverse impact parameter
        double z0 = _params(4);       // longitudinal impact parameter
        
        // Calculate momentum components
        double pT = std::abs(1.0 / qOverPt);
        double theta = 2.0 * std::atan(std::exp(-eta));
        double p = pT / std::sin(theta);
        double charge = (qOverPt > 0) ? particle.charge : -particle.charge;
        
        // 3D momentum vector
        Eigen::Vector3d momentum = this->momentum();
        
        // Calculate beta = v/c for the particle
        double beta = p / std::sqrt(p*p + particle.mass*particle.mass);
        
        // Current position (approximate for the current surface)
        Eigen::Vector3d posAtOrigin = positionAtOrigin();
        
        // Need to calculate position at the current surface
        // This is a simplified approach - in a real implementation,
        // you would need to calculate this more accurately
        dd4hep::rec::Vector3D currentPos = _surface->origin(); // cm to mm
        dd4hep::rec::Vector3D currentDir(momentum(0), momentum(1), momentum(2));
        currentDir = currentDir.unit();
        
        // Calculate intersection with the new surface
        bool intersects = false;
        double pathlength = 0.0;

        // Calculate intersection using surface normal and position
        dd4hep::rec::Vector3D normal = surface()->normal();
        dd4hep::rec::Vector3D surfacePos = surface()->origin();

        // Calculate distance along ray to intersection plane
        double denominator = normal.dot(currentDir);
        if (std::abs(denominator) > 1e-6) { // Avoid divide by zero
            // Calculate distance to plane
            double d = normal.dot(surfacePos - currentPos) / denominator;
            
            // Check if intersection is forward along ray
            if (d >= 0) {
                // Calculate intersection point
                dd4hep::rec::Vector3D intersection(
                    currentPos.x() + currentDir.x() * d,
                    currentPos.y() + currentDir.y() * d,
                    currentPos.z() + currentDir.z() * d);
            
                // Check if point is within surface bounds
                dd4hep::rec::Vector2D localPos = surface()->globalToLocal(intersection);

                // Then check if this position is within bounds
                bool isInside = surface()->insideBounds(intersection);  
            
                if (isInside) {
                    intersects = true;
                    pathlength = d;
                }
            }
        }
        
        if (!intersects) {
            // Handle no intersection case
            return *this;
        }
        
        // Calculate new position
        dd4hep::rec::Vector3D newPos = currentPos + pathlength * currentDir;
        
        // In a real implementation, we would account for the magnetic field
        // by integrating the equations of motion (e.g., using Runge-Kutta)
        // This would give us updated momentum and position at the new surface
        
        // For this simplified version, we'll assume straight line propagation
        // but still update the track parameters for the new position
        
        // Calculate material effects between the two points
        const dd4hep::rec::MaterialVec& materials = 
            const_cast<dd4hep::rec::MaterialManager&>(matMgr).materialsBetween(currentPos, newPos);
        bool isMaterial = !materials.empty();
        double totalThickness = 0.0;
        double radiationLength = 0.0;

        if (isMaterial) {
            // Process all materials along the path
            for (const auto& material_pair : materials) {
                // Each entry is a pair of (Material, thickness)
                const dd4hep::Material& material = material_pair.first;
                double thickness = material_pair.second;
        
                totalThickness += thickness;
        
                // Access radiation length from the Material object directly
                double radLength = material.radLength(); // Or radiationLength() - check which exists
                if (radLength > 0) {
                    radiationLength += thickness / radLength;
                }
            }
        
            // If you need to work with a single effective radiation length
            double effectiveRadLength = (radiationLength > 0) ? totalThickness / radiationLength : 0.0;
        }

        // Initialize transport with identity matrix to avoid numerical issues
        Eigen::Matrix5d transportMatrix = Eigen::Matrix5d::Identity();
        
        // Transport covariance matrix using the transport matrix (jacobian)
        Eigen::Matrix5d transportedCov = transportMatrix * _cov * transportMatrix.transpose();
        
        // Initialize new parameters with the current ones
        Eigen::Vector5d newParams = _params;
        
        // Apply material effects if material is found along the path
        if (isMaterial && radiationLength > 0.0) {
            // Use the effectiveRadLength we calculated earlier
            double effectiveRadLength = (radiationLength > 0) ? totalThickness / radiationLength : 0.0;
        
            // Calculate multiple scattering angle using Highland formula
            double theta0 = 13.6 / (beta * p) * std::sqrt(effectiveRadLength) * (1.0 + 0.038 * std::log(effectiveRadLength));
        
            // Add multiple scattering contribution to covariance matrix
            addMultipleScatteringNoise(transportedCov, theta0);
        
            // Calculate energy loss
            double dE = calculateEnergyLoss(p, beta, effectiveRadLength, particle);
        
            // Update momentum due to energy loss
            double newP = p - dE;
            if (newP > 0) {
                // Recalculate q/pT with the new momentum
                double newPt = newP * std::sin(theta);
                newParams(0) = (qOverPt > 0 ? 1.0 : -1.0) / newPt;
            }
        }

        // Check for NaN or infinities in parameters and covariance
        for (int i = 0; i < 5; ++i) {
            if (std::isnan(newParams(i)) || std::isinf(newParams(i))) {
                // Replace with reasonable values if they're invalid
                newParams(i) = _params(i);  // Keep original parameters
            }
            
            for (int j = 0; j < 5; ++j) {
                if (std::isnan(transportedCov(i,j)) || std::isinf(transportedCov(i,j))) {
                    // Replace with reasonable values if they're invalid
                    transportedCov(i,j) = (i == j) ? 1.0 : 0.0;  // Use identity matrix as fallback
                }
            }
        }
        
        // Return new state
        return TrackState(newParams, transportedCov, newSurface);
    } catch (const std::exception& ex) {
        // Handle the error and log useful diagnostic information
        std::stringstream ss;
        ss << "Error in predictTo: " << ex.what() 
           << " - Surface at (" << newSurface->origin().x() 
           << ", " << newSurface->origin().y() 
           << ", " << newSurface->origin().z() << ")";
        
        // Create a safe fallback state
        // Just return the original state but associated with the new surface
        return TrackState(_params, _cov, newSurface);
    } catch (...) {
        // Catch any other type of exception
        // Create a safe fallback state
        return TrackState(_params, _cov, newSurface);
    }
}


void TrackState::addMultipleScatteringNoise(Eigen::Matrix5d& cov, double theta0) const {
    // Multiple scattering affects primarily the angular parameters (phi, eta)
    
    // Add to the diagonal elements for angular uncertainty
    cov(1, 1) += theta0 * theta0;  // phi variance
    cov(2, 2) += theta0 * theta0 / (std::sinh(_params(2)) * std::sinh(_params(2)));  // eta variance
    
    // Add correlation between phi and eta
    double correlation = theta0 * theta0 * 0.5;  // Simplified correlation factor
    cov(1, 2) += correlation;
    cov(2, 1) += correlation;
    
    // Also add correlations to the impact parameters (d0, z0)
    // This accounts for the displacement due to multiple scattering
    double distance = 1.0;  // Placeholder - should be the distance to the next measurement
    
    cov(1, 3) += theta0 * theta0 * distance;  // phi-d0 correlation
    cov(3, 1) += theta0 * theta0 * distance;
    cov(2, 4) += theta0 * theta0 * distance;  // eta-z0 correlation
    cov(4, 2) += theta0 * theta0 * distance;
}

double TrackState::calculateEnergyLoss(double momentum, double beta, double radLength, 
                                      const ParticleProperties& particle) const {
    // Energy loss calculation depends on the particle type
    double dE = 0.0;
    
    if (particle.name == "electron" || particle.name == "positron") {
        // Electrons/positrons: mainly bremsstrahlung
        // Simplified Bethe-Heitler formula
        double X0 = 0.35;  // Radiation length for silicon in cm (approximate)
        dE = momentum * (1.0 - std::exp(-radLength / X0));
    } else {
        // Heavier particles: mainly ionization
        // Simplified Bethe-Bloch formula
        double K = 0.307075;  // MeV⋅cm²/g
        double Z_A = 0.5;     // Z/A for silicon (approximate)
        double density = 2.33; // g/cm³ for silicon
        double I = 173e-6;    // Mean excitation energy for silicon in GeV
        
        double gamma = 1.0 / std::sqrt(1.0 - beta*beta);
        double mass = particle.mass;
        
        // Maximum energy transfer in a single collision
        double Wmax = 2.0 * electron_mass * beta*beta * gamma*gamma / 
            (1.0 + 2.0 * gamma * electron_mass / mass + std::pow(electron_mass / mass, 2));
        
        // Mean energy loss rate (dE/dx)
        double dEdx = K * Z_A * density * std::pow(particle.charge / beta, 2) *
            (std::log(2.0 * electron_mass * beta*beta * gamma*gamma * Wmax / (I*I)) - 
             2.0 * beta*beta);
        
        // Convert to GeV and multiply by thickness
        dE = dEdx * radLength * 1e-3;  // Convert MeV to GeV
    }
    
    // Ensure energy loss is not greater than the particle's energy
    dE = std::min(dE, momentum - 0.01);  // Leave some minimum momentum
    
    return dE;
}

TrackState TrackState::update(const edm4hep::TrackerHitPlane& hit, 
                             const dd4hep::rec::Surface* surface) const {
    // Get hit position
    const auto& pos = hit.getPosition();
    
    // Convert to local coordinates on surface
    dd4hep::rec::Vector3D hitPosVec(pos[0], pos[1], pos[2]);
    dd4hep::rec::Vector2D localPos = surface->globalToLocal(dd4hep::mm * hitPosVec);
    
    // Create measurement vector
    Eigen::Vector2d meas(localPos[0] / dd4hep::mm, localPos[1] / dd4hep::mm);
    
    // Create measurement covariance matrix
    Eigen::Matrix2d measCov = Eigen::Matrix2d::Zero();
    
    // Get hit covariance - convert from global 3D to local 2D
    // This is a simplified approach - a real implementation would
    // transform the covariance matrix properly
    
    // Get the rotation part of the transformation
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    
    dd4hep::rec::Vector3D normal = surface->normal();
    dd4hep::rec::Vector3D uVec = surface->u();
    dd4hep::rec::Vector3D vVec = surface->v();

    // Normalize these vectors to ensure they're unit vectors
    dd4hep::rec::Vector3D normNormal(
        normal.x() / normal.r(),
        normal.y() / normal.r(),
        normal.z() / normal.r()
    );

    dd4hep::rec::Vector3D normU(
        uVec.x() / uVec.r(),
        uVec.y() / uVec.r(),
        uVec.z() / uVec.r()
    );

    dd4hep::rec::Vector3D normV(
        vVec.x() / vVec.r(),
        vVec.y() / vVec.r(),
        vVec.z() / vVec.r()
    );

    // Fill rotation matrix with basis vectors
    rotation(0,0) = normU.x(); rotation(0,1) = normV.x(); rotation(0,2) = normNormal.x();
    rotation(1,0) = normU.y(); rotation(1,1) = normV.y(); rotation(1,2) = normNormal.y();
    rotation(2,0) = normU.z(); rotation(2,1) = normV.z(); rotation(2,2) = normNormal.z();
    
    // Construct 3D covariance matrix from hit
    Eigen::Matrix3d cov3d = Eigen::Matrix3d::Zero();
    const auto& covValues = hit.getCovMatrix();
    
    // EDM4hep stores covariance as a 6-element vector [xx, xy, xz, yy, yz, zz]
    cov3d(0, 0) = covValues[0]; // xx
    cov3d(0, 1) = covValues[1]; // xy
    cov3d(0, 2) = covValues[2]; // xz
    cov3d(1, 0) = covValues[1]; // xy
    cov3d(1, 1) = covValues[3]; // yy
    cov3d(1, 2) = covValues[4]; // yz
    cov3d(2, 0) = covValues[2]; // xz
    cov3d(2, 1) = covValues[4]; // yz
    cov3d(2, 2) = covValues[5]; // zz
    
    // Project to local coordinates
    dd4hep::rec::Vector3D u_axis = surface->u();
    dd4hep::rec::Vector3D v_axis = surface->v();
    
    Eigen::Matrix<double, 2, 3> projection;
    projection << u_axis.x(), u_axis.y(), u_axis.z(),
                  v_axis.x(), v_axis.y(), v_axis.z();
    
    measCov = projection * cov3d * projection.transpose();
    
    // Project track state to measurement space
    Eigen::Matrix<double, 2, 5> H = projectionMatrix(surface);
    
    // Predicted measurement
    Eigen::Vector2d predictedMeas = H * _params;
    
    // Innovation (residual)
    Eigen::Vector2d residual = meas - predictedMeas;
    
    // Innovation covariance
    Eigen::Matrix2d S = H * _cov * H.transpose() + measCov;
    
    // Kalman gain
    Eigen::Matrix<double, 5, 2> K = _cov * H.transpose() * S.inverse();
    
    // Updated state parameters
    Eigen::Vector5d updatedParams = _params + K * residual;
    
    // Updated covariance
    Eigen::Matrix5d updatedCov = (Eigen::Matrix5d::Identity() - K * H) * _cov;
    
    // Return updated state
    return TrackState(updatedParams, updatedCov, surface);
}

double TrackState::getChi2Increment(const edm4hep::TrackerHitPlane& hit, 
                                   const dd4hep::rec::Surface* surface) const {
    std::cout << "==== Chi2 Calculation Debug ====" << std::endl;
    
    // Get hit position
    const auto& pos = hit.getPosition();
    std::cout << "Hit global position (mm): (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    
    // Convert to local coordinates on surface
    dd4hep::rec::Vector3D hitPosVec(pos[0], pos[1], pos[2]);

    // Try to print surface normal and basis vectors
    dd4hep::rec::Vector3D origin = surface->origin(); //cm to mm
    dd4hep::rec::Vector3D normal = surface->normal();
    dd4hep::rec::Vector3D uVec = surface->u();
    dd4hep::rec::Vector3D vVec = surface->v();

    std::cout << "Surface origin: (" << origin.x() << ", " << origin.y() << ", " << origin.z() << ")" << std::endl;
    std::cout << "Surface normal: (" << normal.x() << ", " << normal.y() << ", " << normal.z() << ")" << std::endl;
    std::cout << "Surface u vector: (" << uVec.x() << ", " << uVec.y() << ", " << uVec.z() << ")" << std::endl;
    std::cout << "Surface v vector: (" << vVec.x() << ", " << vVec.y() << ", " << vVec.z() << ")" << std::endl;
    
    dd4hep::rec::Vector2D localPos = surface->globalToLocal(dd4hep::mm * hitPosVec);
    std::cout << "Hit local position (mm): (" << localPos.u() / dd4hep::mm << ", " 
            << localPos.v() / dd4hep::mm << ")" << std::endl;
    
    // Create measurement vector - local coordinates in mm
    Eigen::Vector2d meas(localPos.u() / dd4hep::mm, localPos.v() / dd4hep::mm);
    std::cout << "Measurement vector: (" << meas(0) << ", " << meas(1) << ")" << std::endl;
    
    // Create measurement covariance matrix
    Eigen::Matrix2d measCov = Eigen::Matrix2d::Zero();
    
    // Get hit covariance - convert from global 3D to local 2D
    // This is a simplified approach - a real implementation would
    // transform the covariance matrix properly
    
    // Get the rotation part of the transformation
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    
    dd4hep::rec::Vector3D normNormal(
        normal.x() / normal.r(),
        normal.y() / normal.r(),
        normal.z() / normal.r()
    );

    dd4hep::rec::Vector3D normU(
        uVec.x() / uVec.r(),
        uVec.y() / uVec.r(),
        uVec.z() / uVec.r()
    );

    dd4hep::rec::Vector3D normV(
        vVec.x() / vVec.r(),
        vVec.y() / vVec.r(),
        vVec.z() / vVec.r()
    );

    // Fill rotation matrix with basis vectors
    rotation(0,0) = normU.x(); rotation(0,1) = normV.x(); rotation(0,2) = normNormal.x();
    rotation(1,0) = normU.y(); rotation(1,1) = normV.y(); rotation(1,2) = normNormal.y();
    rotation(2,0) = normU.z(); rotation(2,1) = normV.z(); rotation(2,2) = normNormal.z();
    std::cout << "Rotation matrix:" << std::endl;
    std::cout << "  [" << rotation(0,0) << ", " << rotation(0,1) << ", " << rotation(0,2) << "]" << std::endl;
    std::cout << "  [" << rotation(1,0) << ", " << rotation(1,1) << ", " << rotation(1,2) << "]" << std::endl;
    std::cout << "  [" << rotation(2,0) << ", " << rotation(2,1) << ", " << rotation(2,2) << "]" << std::endl;
    
    // Construct 3D covariance matrix from hit
    Eigen::Matrix3d cov3d = Eigen::Matrix3d::Zero();
    const auto& covValues = hit.getCovMatrix();
    
    // EDM4hep stores covariance as a 6-element vector [xx, xy, xz, yy, yz, zz]
    cov3d(0, 0) = covValues[0]; // xx
    cov3d(0, 1) = covValues[1]; // xy
    cov3d(0, 2) = covValues[2]; // xz
    cov3d(1, 0) = covValues[1]; // xy
    cov3d(1, 1) = covValues[3]; // yy
    cov3d(1, 2) = covValues[4]; // yz
    cov3d(2, 0) = covValues[2]; // xz
    cov3d(2, 1) = covValues[4]; // yz
    cov3d(2, 2) = covValues[5]; // zz
    std::cout << "3D covariance matrix:" << std::endl;
    std::cout << "  [" << cov3d(0,0) << ", " << cov3d(0,1) << ", " << cov3d(0,2) << "]" << std::endl;
    std::cout << "  [" << cov3d(1,0) << ", " << cov3d(1,1) << ", " << cov3d(1,2) << "]" << std::endl;
    std::cout << "  [" << cov3d(2,0) << ", " << cov3d(2,1) << ", " << cov3d(2,2) << "]" << std::endl;
    
    // Project to local coordinates
    Eigen::Matrix<double, 2, 3> projection;
    projection << uVec.x(), uVec.y(), uVec.z(),
                  vVec.x(), vVec.y(), vVec.z();
    
    measCov = projection * cov3d * projection.transpose();
    std::cout << "Measurement covariance matrix:" << std::endl;
    std::cout << "  [" << measCov(0,0) << ", " << measCov(0,1) << "]" << std::endl;
    std::cout << "  [" << measCov(1,0) << ", " << measCov(1,1) << "]" << std::endl;
    
    // If the covariance matrix has zeros on the diagonal, add a reasonable value
    if (measCov(0, 0) < 1e-6 || measCov(1, 1) < 1e-6) {
        std::cout << "WARNING: Measurement covariance is too small, adding minimum value" << std::endl;
        // Add a reasonable detector resolution - typically ~0.1mm for silicon detectors
        const double minResolution = 0.1 * 0.1;  // 0.1mm resolution squared
        measCov(0, 0) = std::max(measCov(0, 0), minResolution);
        measCov(1, 1) = std::max(measCov(1, 1), minResolution);
        
        std::cout << "Updated measurement covariance:" << std::endl;
        std::cout << "  [" << measCov(0,0) << ", " << measCov(0,1) << "]" << std::endl;
        std::cout << "  [" << measCov(1,0) << ", " << measCov(1,1) << "]" << std::endl;
    }
    
    // Extract track parameters
    double qOverPt = _params(0);
    double phi = _params(1);
    double eta = _params(2);
    double d0 = _params(3);
    double z0 = _params(4);
    
    std::cout << "Track state parameters: (" 
            << qOverPt << ", "  // q/pT
            << phi << ", "  // phi
            << eta << ", "  // eta
            << d0 << ", "  // d0
            << z0 << ")" << std::endl;  // z0
    
    // Calculate theta from eta
    double theta = 2.0 * std::atan(std::exp(-eta));
    
    // Calculate track direction
    Eigen::Vector3d dir;
    dir << std::cos(phi) * std::sin(theta),
           std::sin(phi) * std::sin(theta),
           std::cos(theta);
    
    std::cout << "Track direction: (" << dir(0) << ", " << dir(1) << ", " << dir(2) << ")" << std::endl;
    
    // Calculate track reference position
    Eigen::Vector3d trackPos;
    trackPos << -d0 * std::sin(phi),
                d0 * std::cos(phi),
                z0;
    
    std::cout << "Track reference position: (" << trackPos(0) << ", " << trackPos(1) << ", " << trackPos(2) << ")" << std::endl;
    
    // Calculate the intersection of the track with the surface
    dd4hep::rec::Vector3D surfNormal = surface->normal().unit();
    dd4hep::rec::Vector3D surfOrigin = surface->origin(); // cm to mm
    
    double t = (surfNormal.dot(surfOrigin) - surfNormal.dot(dd4hep::rec::Vector3D(trackPos.x(), trackPos.y(), trackPos.z()))) / 
              surfNormal.dot(dd4hep::rec::Vector3D(dir.x(), dir.y(), dir.z()));
    
    std::cout << "Distance to intersection: " << t << std::endl;
    
    // Check if the intersection is valid
    if (t <= 0) {
        std::cout << "WARNING: No forward intersection, calculating chi2 using projection matrix" << std::endl;
        
        // Use projection matrix method as fallback
        Eigen::Matrix<double, 2, 5> H = projectionMatrix(surface);
        std::cout << "Projection matrix H:" << std::endl;
        std::cout << "  [" << H(0,0) << ", " << H(0,1) << ", " << H(0,2) << ", " << H(0,3) << ", " << H(0,4) << "]" << std::endl;
        std::cout << "  [" << H(1,0) << ", " << H(1,1) << ", " << H(1,2) << ", " << H(1,3) << ", " << H(1,4) << "]" << std::endl;
        
        // Predicted measurement
        Eigen::Vector2d predictedMeas = H * _params;
        std::cout << "Predicted measurement: (" << predictedMeas(0) << ", " << predictedMeas(1) << ")" << std::endl;
        
        // Innovation (residual)
        Eigen::Vector2d residual = meas - predictedMeas;
        std::cout << "Residual (measurement - prediction): (" << residual(0) << ", " << residual(1) << ")" << std::endl;
        
        // Innovation covariance
        Eigen::Matrix2d S = H * _cov * H.transpose() + measCov;
        std::cout << "Innovation covariance matrix S:" << std::endl;
        std::cout << "  [" << S(0,0) << ", " << S(0,1) << "]" << std::endl;
        std::cout << "  [" << S(1,0) << ", " << S(1,1) << "]" << std::endl;
        
        // Check if S is invertible
        Eigen::FullPivLU<Eigen::Matrix2d> lu_decomp(S);
        if (!lu_decomp.isInvertible()) {
            std::cout << "WARNING: Innovation covariance is not invertible!" << std::endl;
            // Add a small diagonal term to make it invertible
            S(0, 0) += 0.1;
            S(1, 1) += 0.1;
        }
        
        // Chi-square contribution
        double chi2 = residual.transpose() * S.inverse() * residual;
        std::cout << "Chi-square value: " << chi2 << std::endl;
        
        std::cout << "==== End Chi2 Calculation Debug ====" << std::endl;
        return chi2;
    }
    
    // Calculate intersection point
    Eigen::Vector3d intersectPos = trackPos + t * dir;
    std::cout << "Intersection point: (" << intersectPos(0) << ", " << intersectPos(1) << ", " << intersectPos(2) << ")" << std::endl;
    
    // Convert to local coordinates
    dd4hep::rec::Vector3D globalPos(intersectPos.x(), intersectPos.y(), intersectPos.z());
    dd4hep::rec::Vector2D predLocalPos = surface->globalToLocal(globalPos);
    
    // Create predicted measurement vector
    Eigen::Vector2d predictedMeas(predLocalPos.u() / dd4hep::mm, predLocalPos.v() / dd4hep::mm);
    std::cout << "Predicted measurement: (" << predictedMeas(0) << ", " << predictedMeas(1) << ")" << std::endl;
    
    // Calculate residual
    Eigen::Vector2d residual = meas - predictedMeas;
    std::cout << "Residual (measurement - prediction): (" << residual(0) << ", " << residual(1) << ")" << std::endl;
    
    // Calculate the innovation covariance
    // Still use the projection matrix to relate parameter uncertainties to measurement uncertainties
    Eigen::Matrix<double, 2, 5> H = projectionMatrix(surface);
    Eigen::Matrix2d paramCov = H * _cov * H.transpose();
    
    std::cout << "Parameter covariance matrix:" << std::endl;
    std::cout << "  [" << paramCov(0,0) << ", " << paramCov(0,1) << "]" << std::endl;
    std::cout << "  [" << paramCov(1,0) << ", " << paramCov(1,1) << "]" << std::endl;
    
    // Add a minimum uncertainty from track prediction
    if (paramCov(0, 0) < 1e-6 || paramCov(1, 1) < 1e-6) {
        std::cout << "WARNING: Parameter covariance is too small, adding minimum value" << std::endl;
        const double minTrackUncertainty = 1.0; // 1mm uncertainty in track prediction
        paramCov(0, 0) = std::max(paramCov(0, 0), minTrackUncertainty);
        paramCov(1, 1) = std::max(paramCov(1, 1), minTrackUncertainty);
        
        std::cout << "Updated parameter covariance:" << std::endl;
        std::cout << "  [" << paramCov(0,0) << ", " << paramCov(0,1) << "]" << std::endl;
        std::cout << "  [" << paramCov(1,0) << ", " << paramCov(1,1) << "]" << std::endl;
    }
    
    // Innovation covariance
    Eigen::Matrix2d S = paramCov + measCov;
    std::cout << "Innovation covariance matrix S:" << std::endl;
    std::cout << "  [" << S(0,0) << ", " << S(0,1) << "]" << std::endl;
    std::cout << "  [" << S(1,0) << ", " << S(1,1) << "]" << std::endl;
    
    // Ensure matrix is positive definite
    Eigen::FullPivLU<Eigen::Matrix2d> lu_decomp(S);
    if (!lu_decomp.isInvertible()) {
        std::cout << "WARNING: Innovation covariance is not invertible!" << std::endl;
        // Add a small diagonal term to make it invertible
        S(0, 0) += 0.1;
        S(1, 1) += 0.1;
    }
    
    // Chi-square contribution
    double chi2 = residual.transpose() * S.inverse() * residual;
    
    // Calculate individual terms to see where the large values come from
    // We can decompose chi2 = r^T * S^-1 * r into components
    Eigen::Matrix2d Sinv = S.inverse();
    std::cout << "S inverse:" << std::endl;
    std::cout << "  [" << Sinv(0,0) << ", " << Sinv(0,1) << "]" << std::endl;
    std::cout << "  [" << Sinv(1,0) << ", " << Sinv(1,1) << "]" << std::endl;
    
    double term1 = residual(0) * Sinv(0,0) * residual(0);
    double term2 = residual(0) * Sinv(0,1) * residual(1);
    double term3 = residual(1) * Sinv(1,0) * residual(0);
    double term4 = residual(1) * Sinv(1,1) * residual(1);
    
    std::cout << "Chi2 components:" << std::endl;
    std::cout << "  Term 1 (residual[0] * Sinv[0,0] * residual[0]): " << term1 << std::endl;
    std::cout << "  Term 2 (residual[0] * Sinv[0,1] * residual[1]): " << term2 << std::endl;
    std::cout << "  Term 3 (residual[1] * Sinv[1,0] * residual[0]): " << term3 << std::endl;
    std::cout << "  Term 4 (residual[1] * Sinv[1,1] * residual[1]): " << term4 << std::endl;
    std::cout << "  Sum: " << term1 + term2 + term3 + term4 << " (should equal chi2)" << std::endl;
    
    std::cout << "Chi-square value: " << chi2 << std::endl;
    
    std::cout << "==== End Chi2 Calculation Debug ====" << std::endl;
    
    return chi2;
}

//------------------------------------------------------------------------------
// Track Implementation
//------------------------------------------------------------------------------

Track::Track(const TrackState& initialState) 
    : _chi2(0.0) {
    _states.push_back(initialState);
}

void Track::addHit(const edm4hep::TrackerHitPlane& hit, 
                  const dd4hep::rec::Surface* surface,
                  const TrackState& state) {
    _hits.push_back(hit);
    _surfaces.push_back(surface);
    _states.push_back(state);
    
    // Update chi-square
    _chi2 += state.getChi2Increment(hit, surface);
}

Eigen::Vector3d Track::momentumAt(const dd4hep::Position& point) const {
    // Find the closest state to the given point
    size_t closestIdx = 0;
    double minDist = std::numeric_limits<double>::max();
    
    for (size_t i = 0; i < _states.size(); ++i) {
        // Get state position
        Eigen::Vector3d statePos = _states[i].positionAtOrigin();
        
        double dist = (point.x() - statePos.x()) * (point.x() - statePos.x()) + 
                      (point.y() - statePos.y()) * (point.y() - statePos.y()) + 
                      (point.z() - statePos.z()) * (point.z() - statePos.z());
        
        if (dist < minDist) {
            minDist = dist;
            closestIdx = i;
        }
    }
    
    // Get momentum from the closest state
    return _states[closestIdx].momentum();
}

//------------------------------------------------------------------------------
// KalmanTracking Implementation
//------------------------------------------------------------------------------

// Constructor
KalmanTracking::KalmanTracking(const std::string& name, ISvcLocator* pSvcLocator)
    : Algorithm(name, pSvcLocator),
      m_particleProperties(muon_mass, -1.0, "muon") {
    
    // Initialize particle map with common particles
    m_particleMap["electron"] = ParticleProperties(electron_mass, -1.0, "electron");
    m_particleMap["positron"] = ParticleProperties(electron_mass, 1.0, "positron");
    m_particleMap["muon"] = ParticleProperties(muon_mass, -1.0, "muon");
    m_particleMap["antimuon"] = ParticleProperties(muon_mass, 1.0, "antimuon");
    m_particleMap["pion"] = ParticleProperties(pion_mass, 1.0, "pion");
    m_particleMap["kaon"] = ParticleProperties(kaon_mass, 1.0, "kaon");
    m_particleMap["proton"] = ParticleProperties(proton_mass, 1.0, "proton");
    m_particleMap["antiproton"] = ParticleProperties(proton_mass, -1.0, "antiproton");
}

// Initialize method
StatusCode KalmanTracking::initialize() {
    StatusCode sc = Algorithm::initialize();
    if (!sc.isSuccess()) return sc;

    // Get detector service
    if (!m_geoSvc.retrieve()) {
        error() << "Failed to retrieve GeoSvc" << endmsg;
        return StatusCode::FAILURE;
    }

    // Get detector
    m_detector = m_geoSvc->getDetector();
    if (!m_detector) {
        error() << "Failed to access detector" << endmsg;
        return StatusCode::FAILURE;
    }

    // Get BitFieldCoder for cell ID decoding
    std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringParameter);
    m_bitFieldCoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDEncodingString);
    
    // Get surface map
    auto surfaceMan = m_detector->extension<dd4hep::rec::SurfaceManager>();
    m_surfaceMap = surfaceMan->map(m_detectorName);
    if (!m_surfaceMap) {
        error() << "Could not find surface map for detector: " << m_detectorName << endmsg;
        return StatusCode::FAILURE;
    }

    // Get magnetic field
    m_field = m_detector->field();
    
    // Create material manager
    dd4hep::Volume worldVol = m_detector->worldVolume();
    m_materialManager = new dd4hep::rec::MaterialManager(worldVol);
    //m_materialManager = new dd4hep::rec::MaterialManager(*m_detector);

    // Get all surfaces and group them by layer 
    // Another possible and simpler solution is : to call surface by layer directly from surfaceMap,
    // but for the moment we keep the following effeicient way of surfaceList.
    dd4hep::DetElement world = m_detector->world();
    dd4hep::rec::SurfaceHelper surfHelper(world);
    //dd4hep::rec::SurfaceHelper surfHelper(*m_detector);

    const dd4hep::rec::SurfaceList& surfList = surfHelper.surfaceList();

    // Convert to your vector of Surface pointers
    m_surfaces.clear();
    for (const auto& surf : surfList) {
        // Cast to Surface* if necessary
        auto surface = dynamic_cast<const dd4hep::rec::Surface*>(surf);
        if (surface) {
            m_surfaces.push_back(surface);
        }
    }
    
    // Group surfaces by layer
    m_surfacesByLayer = getSurfacesByLayer();
    
    // Set particle properties based on configuration
    auto it = m_particleMap.find(m_particleType);
    if (it != m_particleMap.end()) {
        m_particleProperties = it->second;
        info() << "Using particle type: " << m_particleType 
               << " (mass=" << m_particleProperties.mass 
               << " GeV, charge=" << m_particleProperties.charge << ")" << endmsg;
    } else {
        warning() << "Unknown particle type: " << m_particleType 
                 << ", defaulting to pion" << endmsg;
    }
    
    info() << "Found " << m_surfaces.size() << " detector surfaces in " 
           << m_surfacesByLayer.size() << " layers" << endmsg;

    return StatusCode::SUCCESS;
}

StatusCode KalmanTracking::configure() {
    // No additional configuration needed beyond what's in initialize()
    return StatusCode::SUCCESS;
}

// Execute method
StatusCode KalmanTracking::execute() {
    // Retrieve input hit collection
    DataObject* obj = nullptr;
    StatusCode sc = eventSvc()->retrieveObject(m_inputHitCollection.value(), obj);

    if (sc.isFailure() || obj == nullptr) {
        error() << "Failed to retrieve collection: " << m_inputHitCollection.value() << endmsg;
        return sc;
    }

    // Try to get the collection base from the wrapper
    auto* wrapper = dynamic_cast<AnyDataWrapper<std::unique_ptr<podio::CollectionBase>>*>(obj);
    if (!wrapper) {
        error() << "Retrieved object is not a podio collection wrapper" << endmsg;
        return StatusCode::FAILURE;
    }

    // Get the collection base
    auto& colBasePtr = wrapper->getData();
    podio::CollectionBase* colBase = colBasePtr.get();

    if (!colBase) {
        error() << "Null collection base pointer" << endmsg;
        return StatusCode::FAILURE;
    }

    // Now try to cast to the specific collection type
    edm4hep::TrackerHitPlaneCollection* hitCollection = 
        static_cast<edm4hep::TrackerHitPlaneCollection*>(colBase);

    if (!hitCollection) {
        error() << "Failed to cast input collection to TrackerHitPlaneCollection" << endmsg;
        return StatusCode::FAILURE;
    }

    // Try to access the event header to get event number
    DataObject* evtObj = nullptr;
    sc = eventSvc()->retrieveObject("EventHeader", evtObj);
    if (sc.isSuccess() && evtObj != nullptr) {
        // Try to get the collection base from the wrapper
        auto* wrapper = dynamic_cast<AnyDataWrapper<std::unique_ptr<podio::CollectionBase>>*>(evtObj);
        if (wrapper) {
            auto& colBasePtr = wrapper->getData();
            podio::CollectionBase* colBase = colBasePtr.get();
            if (colBase) {
                // Now try to cast to the specific collection type
                auto* evtHeader = static_cast<edm4hep::EventHeaderCollection*>(colBase);
                if (evtHeader && evtHeader->size() > 0) {
                    debug() << "  Event number: " << (*evtHeader)[0].getEventNumber() << endmsg;
                } else {
                    debug() << "  Could not get event number from collection" << endmsg;
                }
            } else {
                debug() << "  Null collection base pointer for EventHeader" << endmsg;
            }
        } else {
            debug() << "  Retrieved EventHeader is not a podio collection wrapper" << endmsg;
        }
    } else {
        debug() << "  Could not get event number" << endmsg;
    }

    info() << "Processing " << hitCollection->size() << " tracker hits" << endmsg;
    // getting info about first hits 
    if (hitCollection->size() > 0) {
        debug() << "First hit details:" << endmsg;
        const auto& hit = (*hitCollection)[0];
        debug() << "  CellID: " << hit.getCellID() << endmsg;
        debug() << "  Position: (" << hit.getPosition()[0] << ", " 
                                  << hit.getPosition()[1] << ", " 
                                  << hit.getPosition()[2] << ")" << endmsg;
    } else {
        debug() << "WARNING: Negative number of hits!!!" << endmsg;
    }
    
    // Check if we have enough hits for tracking
    std::vector<Track> tracks;
    if (hitCollection->size() < 3) {
        warning() << "Not enough hits to create tracks. Need at least 3." << endmsg;
        // Create empty track vector
    } else {
        // Only try to find tracks if we have enough hits
        try {
            tracks = findTracks(hitCollection);
        } catch (const std::exception& ex) {
            error() << "Exception during track finding: " << ex.what() << endmsg;
            // Continue with empty track vector
        }
    }

    info() << "Found " << tracks.size() << " tracks" << endmsg;
    
    // Create output track collection with proper initialization
    auto* trackCollection = new edm4hep::TrackCollection();

    // Make sure to register using the correct approach for EDM4hep
    // Often with EDM4hep, collections need special handling
    DataObject* trackDataObject = nullptr;

    // Try different approaches to convert/wrap the collection
    try {
        // Option 1: Direct cast (might not work)
        trackDataObject = dynamic_cast<DataObject*>(trackCollection);
        
        if (!trackDataObject) {
            // Option 2: Create a wrapper (common in EDM4hep/podio)
            trackDataObject = new AnyDataWrapper<std::unique_ptr<podio::CollectionBase>>(
                std::unique_ptr<podio::CollectionBase>(trackCollection));
        }
    } catch (const std::exception& ex) {
        error() << "Exception creating output collection: " << ex.what() << endmsg;
        delete trackCollection;
        return StatusCode::FAILURE;
    }

    if (!trackDataObject) {
        error() << "Failed to convert/wrap TrackCollection to DataObject" << endmsg;
        delete trackCollection;
        return StatusCode::FAILURE;
    }

    // Convert internal track representation to standard track format
    for (const auto& track : tracks) {
        createTrack(trackCollection, track);
    }
    
    // Register output collection
   //DataObject* trackDataObject = dynamic_cast<DataObject*>(trackCollection);
    if (!trackDataObject) {
        error() << "Failed to convert TrackCollection to DataObject" << endmsg;
        delete trackCollection;
        return StatusCode::FAILURE;
    }
    sc = eventSvc()->registerObject(m_outputTrackCollection.value(), trackDataObject);
    if (sc.isFailure()) {
        error() << "Failed to register output track collection" << endmsg;
        delete trackCollection;  // Prevent memory leak
        return sc;
    }
    
    return StatusCode::SUCCESS;
}

// Finalize method
StatusCode KalmanTracking::finalize() {
    // Clean up material manager
    if (m_materialManager) {
        delete m_materialManager;
        m_materialManager = nullptr;
    }
    
    return Algorithm::finalize();
}

// Find surface for a hit
const dd4hep::rec::Surface* KalmanTracking::findSurface(const edm4hep::TrackerHitPlane& hit) const {  // should it be edm4hep::TrackerHitPlane!?
    return findSurfaceByID(hit.getCellID());
}

// Find surface by cell ID
const dd4hep::rec::Surface* KalmanTracking::findSurfaceByID(uint64_t cellID) const {
    // Find surface for this cell ID using the surface map
    dd4hep::rec::SurfaceMap::const_iterator sI = m_surfaceMap->find(cellID);
    
    if (sI != m_surfaceMap->end()) {
        return dynamic_cast<const dd4hep::rec::Surface*>(sI->second);
    }
    
    return nullptr;
}

// Extract layer ID from cell ID
int KalmanTracking::getLayerID(uint64_t cellID) const {
    // Use BitFieldCoder to extract layer ID
    return m_bitFieldCoder->get(cellID, "layer");
}

// Group surfaces by detector layer
std::map<int, std::vector<const dd4hep::rec::Surface*>> KalmanTracking::getSurfacesByLayer() const {
    std::map<int, std::vector<const dd4hep::rec::Surface*>> result;
    
    // Process all surfaces and group them by layer ID
    for (const auto& surface : m_surfaces) {
        // Get detector element
        dd4hep::DetElement det = surface->detElement(); // not sure if it detElement ot  detector();
        if (!det.isValid()) continue;
        
        // Try to get layer ID from volume ID
        uint64_t volID = det.volumeID();
        int layerID = 0;
        
        try {
            layerID = m_bitFieldCoder->get(volID, "layer");
        } catch (...) {
            // If no layer field, try to determine from position
            dd4hep::rec::Vector3D pos = surface->origin(); // cm to mm
            layerID = static_cast<int>(std::sqrt(pos.x()*pos.x() + pos.y()*pos.y()) / 10.0);
        }
        
        // Add surface to the appropriate layer group
        result[layerID].push_back(surface);
    }
    
    return result;
}

double KalmanTracking::getRadiationLength(const dd4hep::rec::Vector3D& start, 
                                        const dd4hep::rec::Vector3D& end) const {
    if (!m_materialManager) {
        warning() << "Material manager not initialized for radiation length calculation" << endmsg;
        return 0.0;
    }
    
    // Get materials along the path
    const dd4hep::rec::MaterialVec& materials = 
        m_materialManager->materialsBetween(start, end);
    
    double totalRadiationLength = 0.0;
    
    // Process all materials along the path
    for (const auto& material_pair : materials) {
        // Each entry is a pair of (Material, thickness)
        double thickness = material_pair.second;
        const dd4hep::Material& material = material_pair.first;
        
        // Get radiation length from the Material object
        double radLength = material.radLength(); // Or use material.radiationLength()
        
        if (radLength > 0) {
            totalRadiationLength += thickness / radLength;
        }
    }
    
    return totalRadiationLength;
}

// Find potentially intersecting surfaces for a track
std::vector<const dd4hep::rec::Surface*> KalmanTracking::findIntersectingSurfaces(
    const TrackState& state, double maxDistance) const {
    debug() << "Finding intersecting surfaces within " << maxDistance << " mm" << endmsg;
    
    std::vector<const dd4hep::rec::Surface*> intersectingSurfaces;
    
    // Get track parameters
    Eigen::Vector3d momentum = state.momentum();
    Eigen::Vector3d pos = state.positionAtOrigin();
    
    debug() << "Track position: (" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")" << endmsg;
    debug() << "Track momentum: (" << momentum.x() << ", " << momentum.y() << ", " << momentum.z() << ")" << endmsg;

    // Normalize direction vector
    Eigen::Vector3d dir = momentum.normalized();

    debug() << "Track direction: (" << dir.x() << ", " << dir.y() << ", " << dir.z() << ")" << endmsg;
    debug() << "Checking " << m_surfacesByLayer.size() << " layers for intersections" << endmsg;
    
    // For each layer, find potential intersections
    for (const auto& [layerID, surfaces] : m_surfacesByLayer) {
        bool layerIntersected = false;
        
        for (const auto& surface : surfaces) {
            // Skip current surface
            if (surface == state.surface()) continue;
            
            // Calculate approximate distance to surface
            dd4hep::rec::Vector3D surfPos = surface->origin(); //cm to mm
            Eigen::Vector3d surfCenter(surfPos.x(), surfPos.y(), surfPos.z());
            
            // Project vector from current position to surface center onto track direction
            Eigen::Vector3d toSurface = surfCenter - pos;
            double projection = toSurface.dot(dir);
            
            // Only consider surfaces ahead in track direction
            if (projection > 0) {
                // Calculate closest approach of track to surface center
                Eigen::Vector3d closestPoint = pos + projection * dir;
                double distance = (closestPoint - surfCenter).norm();
                
                // Check if within distance threshold
                if (distance < maxDistance) {
                    // Check for actual intersection
                    dd4hep::rec::Vector3D trackPos(pos.x(), pos.y(), pos.z());
                    dd4hep::rec::Vector3D trackDir(dir.x(), dir.y(), dir.z());
                    
                    bool intersects = false;
                    double pathlength = 0.0;

                    // Calculate intersection using surface normal and position
                    dd4hep::rec::Vector3D normal = surface->normal();
                    dd4hep::rec::Vector3D surfacePos = surface->origin();

                    // Calculate distance along ray to intersection plane
                    double denominator = normal.dot(trackDir);
                    if (std::abs(denominator) > 1e-6) { // Avoid divide by zero
                    // Calculate distance to plane
                    double d = normal.dot(surfacePos - trackPos) / denominator;
                        
                        // Check if intersection is forward along ray
                        if (d >= 0) {
                            // Calculate intersection point
                            dd4hep::rec::Vector3D intersection(
                                trackPos.x() + trackDir.x() * d,
                                trackPos.y() + trackDir.y() * d,
                                trackPos.z() + trackDir.z() * d);
                            
                            // Check if point is within surface bounds
                            dd4hep::rec::Vector2D localPos = surface->globalToLocal(intersection);

                            // Then check if this position is within bounds
                            bool isInside = surface->insideBounds(intersection);  

                            if (isInside) {
                            intersects = true;
                            pathlength = d;
                            }
                        }
                    }
    
                    if (intersects && pathlength > 0) {
                        intersectingSurfaces.push_back(surface);
                        layerIntersected = true;
                    }
                }
            }
        }
        
        // If we found an intersection with this layer, move to the next layer
        // This assumes layers are ordered and we expect one hit per layer
        if (layerIntersected) {
            continue;
        }
    }
    debug() << "Found " << intersectingSurfaces.size() << " intersecting surfaces" << endmsg;

    return intersectingSurfaces;
}

// Core tracking methods
//------------------------------------------------------------------------------

TrackState KalmanTracking::createSeedState(const edm4hep::TrackerHitPlane& hit1, // not sure if it should be TrackerHitPlane
                                         const dd4hep::rec::Surface* surface1,
                                         const edm4hep::TrackerHitPlane& hit2,  // not sure if it should be TrackerHitPlane
                                         const dd4hep::rec::Surface* surface2) {
    debug() << "Inside createSeedState method" << endmsg;

    // Get hit positions
    const auto& pos1 = hit1.getPosition();
    const auto& pos2 = hit2.getPosition();
    
    // Direction vector
    double dx = pos2[0] - pos1[0];
    double dy = pos2[1] - pos1[1];
    double dz = pos2[2] - pos1[2];

    debug() << "Direction vector: (" << dx << ", " << dy << ", " << dz << ")" << endmsg;
    
    // Calculate track parameters
    
    // Azimuthal angle
    double phi = std::atan2(dy, dx);
    
    // Polar angle
    double theta = std::atan2(std::sqrt(dx*dx + dy*dy), dz);
    
    // Pseudorapidity
    double eta = -std::log(std::tan(theta/2.0));

    debug() << "Calculated angles: phi=" << phi << ", theta=" << theta << ", eta=" << eta << endmsg;
    
    // Impact parameters (simplified - assumes origin as reference)
    // In a real implementation, you'd calculate these more accurately
    double d0 = 0.0;  // Approximation
    double z0 = 0.0;  // Approximation
    
    // Use provided initial momentum estimate for the track
    double pT = m_initialMomentum * std::sin(theta);
    double qOverPt = 1.0 / pT;  // Assume positive charge initially

    debug() << "Track parameters: qOverPt=" << qOverPt << ", d0=" << d0 << ", z0=" << z0 << endmsg;
    
    // Initialize track parameters
    Eigen::Vector5d params;
    params << qOverPt, phi, eta, d0, z0;
    
    // Initial covariance matrix - large uncertainties on all parameters
    Eigen::Matrix5d cov = Eigen::Matrix5d::Identity();
    cov(0, 0) = 1.0;    // q/pT uncertainty
    cov(1, 1) = 0.01;   // phi uncertainty [rad]
    cov(2, 2) = 0.01;   // eta uncertainty
    cov(3, 3) = 10.0;   // d0 uncertainty [mm]
    cov(4, 4) = 10.0;   // z0 uncertainty [mm]

    debug() << "Created initial track state successfully" << endmsg;
    
    return TrackState(params, cov, surface1);
}

bool KalmanTracking::extendTrackCandidate(const TrackState& seedState, 
                                       const edm4hep::TrackerHitPlaneCollection* hits,
                                       std::vector<bool>& usedHits,
                                       std::vector<edm4hep::TrackerHitPlane>& trackHits,
                                       std::vector<const dd4hep::rec::Surface*>& trackSurfaces,
                                       std::vector<size_t>& trackHitIndices) {
    debug() << "Extending track candidate from seed state" << endmsg;
    
    // Get current track state
    TrackState currentState = seedState;
    
    // Find potential surfaces that the track might intersect
    debug() << "Finding potential intersecting surfaces..." << endmsg;
    std::vector<const dd4hep::rec::Surface*> potentialSurfaces = 
        findIntersectingSurfaces(currentState, m_maxDistanceToSurface);
    debug() << "Found " << potentialSurfaces.size() << " potentially intersecting surfaces" << endmsg;
    
    // For each potential surface
    for (const auto& surface : potentialSurfaces) {
        debug() << "  Checking surface at (" << surface->origin().x() << ", "
                << surface->origin().y() << ", " << surface->origin().z() << ")" << endmsg;
        // Skip surfaces already used by this track
        bool surfaceUsed = false;
        for (const auto& surf : trackSurfaces) {
            if (surf == surface) {
                surfaceUsed = true;
                break;
            }
        }
        if (surfaceUsed) continue;
        
        try {
            // Predict state to this surface
            debug() << "  Attempting to predict track state to this surface..." << endmsg;
            TrackState predictedState = currentState.predictTo(surface, m_field, *m_materialManager, m_particleProperties);
            debug() << "  Prediction successful!" << endmsg;
            
            debug() << "Looking for compatible hits on surface at (" 
                    << surface->origin().x() << ", " 
                    << surface->origin().y() << ", " 
                    << surface->origin().z() << ")" << endmsg;


            // Find compatible hits on this surface
            auto compatibleHits = findCompatibleHits(predictedState, hits, usedHits);
            debug() << "  Found " << compatibleHits.size() << " compatible hits" << endmsg;
            
            if (!compatibleHits.empty()) {
                // Use the best compatible hit (lowest chi-square)
                size_t hitIdx = std::get<0>(compatibleHits[0]);
                edm4hep::TrackerHitPlane hit = std::get<1>(compatibleHits[0]);
                const dd4hep::rec::Surface* hitSurface = std::get<2>(compatibleHits[0]);
                
                // Update state with this hit
                TrackState updatedState = predictedState.update(hit, hitSurface);
                
                // Add hit to track
                trackHits.push_back(hit);
                trackSurfaces.push_back(hitSurface);
                trackHitIndices.push_back(hitIdx);
                
                // Update current state
                currentState = updatedState;
            }
        } catch (const std::exception& ex) {
            error() << "  Exception when checking surface: " << ex.what() << endmsg;
            // Continue to next surface instead of crashing
            continue;
        }
    }
    // Track is valid if it has at least 3 hits
    bool trackValid = trackHits.size() >= 3;
    debug() << "Track candidate extension finished with " << trackHits.size() 
            << " hits, valid = " << (trackValid ? "true" : "false") << endmsg;
    
    // Track is valid if it has at least 3 hits
    return trackHits.size() >= 3;
}

std::vector<Track> KalmanTracking::findTracks(const edm4hep::TrackerHitPlaneCollection* hits) {
    std::vector<Track> tracks;
    // Check if we have enough hits for tracking (need at least 3)
    if (hits->size() < 3) {
        info() << "Not enough hits for tracking. Need at least 3, found " << hits->size() << endmsg;
        return tracks; // Return empty track vector
    }
    
    std::vector<bool> usedHits(hits->size(), false);
    
    // Group hits by surface
    std::map<const dd4hep::rec::Surface*, std::vector<std::pair<size_t, edm4hep::TrackerHitPlane>>> hitsBySurface; // should it be TrackerHitPlane!?
    
    // Find surfaces for all hits
    for (size_t i = 0; i < hits->size(); ++i) {
        const auto& hit = (*hits)[i];
        const dd4hep::rec::Surface* surface = findSurface(hit);
        
        if (surface) {
            hitsBySurface[surface].push_back(std::make_pair(i, hit));
        } else {
            debug() << "Could not find surface for hit with cellID " << hit.getCellID() << endmsg;
        }
    }

    info() << "Found " << hitsBySurface.size() << " surfaces with hits" << endmsg;
    if (hitsBySurface.empty()) {
        info() << "No surfaces found for any hits! Check detector configuration." << endmsg;
        return tracks; // Return empty tracks if no surfaces found
    }

    // Print how many hits are on each surface
    debug() << "Hit distribution by surface:" << endmsg;
    for (const auto& [surface, hits] : hitsBySurface) {
        debug() << "  Surface at (" << surface->origin().x() << ", " 
                << surface->origin().y() << ", " << surface->origin().z() 
                << ") has " << hits.size() << " hits" << endmsg;
    }

    info() << "Found surfaces for " << hitsBySurface.size() << " different surfaces" << endmsg;
    // Track seeding strategy:
    // 1. Start with pairs of hits from different layers
    // 2. Build track seeds and extend them
    
    // For each layer pair
    for (auto it1 = m_surfacesByLayer.begin(); it1 != m_surfacesByLayer.end(); ++it1) {
        int layer1 = it1->first;
        const auto& surfaces1 = it1->second;
        
        for (auto it2 = std::next(it1); it2 != m_surfacesByLayer.end(); ++it2) {
            int layer2 = it2->first;
            const auto& surfaces2 = it2->second;
            
            // For each surface in first layer
            for (const auto& surface1 : surfaces1) {
                const auto& hitList1 = hitsBySurface[surface1];
                
                // For each surface in second layer
                for (const auto& surface2 : surfaces2) {
                    const auto& hitList2 = hitsBySurface[surface2];
                    
                    // For each hit pair
                    for (const auto& [idx1, hit1] : hitList1) {
                        if (usedHits[idx1]) continue;
                        
                        for (const auto& [idx2, hit2] : hitList2) {
                            if (usedHits[idx2]) continue;

                            // Inside the nested for loops where you create seed states
                            debug() << "Attempting to create seed from hits on layers " 
                                    << layer1 << " and " << layer2 << endmsg;
                            debug() << "  Surface 1 origin: (" << surface1->origin().x() << ", "
                                    << surface1->origin().y() << ", " << surface1->origin().z() << ")" << endmsg;
                            debug() << "  Surface 2 origin: (" << surface2->origin().x() << ", "
                                    << surface2->origin().y() << ", " << surface2->origin().z() << ")" << endmsg;
                            debug() << "  Hit 1 position: (" << hit1.getPosition()[0] << ", " 
                                    << hit1.getPosition()[1] << ", " << hit1.getPosition()[2] << ")" << endmsg;
                            debug() << "  Hit 2 position: (" << hit2.getPosition()[0] << ", " 
                                    << hit2.getPosition()[1] << ", " << hit2.getPosition()[2] << ")" << endmsg;

                            // Right before creating the seed state
                            debug() << "Creating seed state from these two hits..." << endmsg;

                            try {
                                // Create seed state from these two hits
                                TrackState seedState = createSeedState(hit1, surface1, hit2, surface2);
                                debug() << "Seed state created successfully with parameters: ("
                                        << seedState.parameters()(0) << ", " // q/pT
                                        << seedState.parameters()(1) << ", " // phi
                                        << seedState.parameters()(2) << ", " // eta
                                        << seedState.parameters()(3) << ", " // d0
                                        << seedState.parameters()(4) << ")" << endmsg; // z0

                                // Try to extend the track code...
                            } catch (const std::exception& ex) {
                                error() << "Exception during seed creation: " << ex.what() << endmsg;
                                continue; // Try next hit pair
                            }
                            
                            // Create seed state from these two hits
                            TrackState seedState = createSeedState(hit1, surface1, hit2, surface2);
                            
                            // Find more hits and build track candidate
                            std::vector<edm4hep::TrackerHitPlane> trackHits = {hit1, hit2};
                            std::vector<const dd4hep::rec::Surface*> trackSurfaces = {surface1, surface2};
                            std::vector<size_t> trackHitIndices = {idx1, idx2};
                            
                            // Try to extend the track by adding compatible hits
                            bool trackValid = extendTrackCandidate(seedState, hits, usedHits, trackHits, trackSurfaces, trackHitIndices);
                            
                            if (trackValid) {
                                // Fit final track
                                Track track = fitTrack(trackHits, trackSurfaces, seedState);
                                
                                // Quality cuts
                                if (track.chi2() / track.ndf() < m_maxChi2) {
                                    tracks.push_back(track);
                                    
                                    // Mark hits as used
                                    for (size_t idx : trackHitIndices) {
                                        usedHits[idx] = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return tracks;
}

Track KalmanTracking::fitTrack(const std::vector<edm4hep::TrackerHitPlane>& hits,
                             const std::vector<const dd4hep::rec::Surface*>& surfaces,
                             const TrackState& seedState) {
    // Initialize track with seed state
    Track track(seedState);
    
    // Determine the hit ordering based on distance along track
    
    // Get seed direction
    Eigen::Vector3d seedDir = seedState.momentum().normalized();
    Eigen::Vector3d seedPos = seedState.positionAtOrigin();
    
    // Calculate projection of each hit along track direction
    std::vector<std::tuple<double, edm4hep::TrackerHitPlane, const dd4hep::rec::Surface*>> hitProjections;
    for (size_t i = 0; i < hits.size(); ++i) {
        const auto& pos = hits[i].getPosition();
        Eigen::Vector3d posVec(pos[0], pos[1], pos[2]);
        
        // Project hit position onto track direction
        double projection = (posVec - seedPos).dot(seedDir);
        hitProjections.push_back(std::make_tuple(projection, hits[i], surfaces[i]));
    }
    
    // Sort hits by projection value
    std::sort(hitProjections.begin(), hitProjections.end(), 
              [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
    
    // Forward filter pass (Kalman filter)
    TrackState currentState = seedState;
    
    for (const auto& [projection, hit, surface] : hitProjections) {
        // Predict to the hit surface
        TrackState predictedState = currentState.predictTo(surface, m_field, *m_materialManager, m_particleProperties);
        
        // Update with the hit
        TrackState updatedState = predictedState.update(hit, surface);
        
        // Add hit to track
        track.addHit(hit, surface, updatedState);
        
        // Update current state
        currentState = updatedState;
    }
    
    return track;
}

std::vector<std::tuple<size_t, edm4hep::TrackerHitPlane, const dd4hep::rec::Surface*>> KalmanTracking::findCompatibleHits(
    const TrackState& state, 
    const edm4hep::TrackerHitPlaneCollection* hits,
    const std::vector<bool>& usedHits) {
    
    std::vector<std::tuple<size_t, edm4hep::TrackerHitPlane, const dd4hep::rec::Surface*>> compatibleHits;
    
    // Check each hit for compatibility
    for (size_t i = 0; i < hits->size(); ++i) {
        // Skip already used hits
        if (usedHits[i]) continue;
        
        const auto& hit = (*hits)[i];
        const dd4hep::rec::Surface* surface = findSurface(hit);
        
        // Skip hits with no surface or not on the current surface
        //if (!surface || surface != state.surface()) continue;
        
        // Calculate chi-square
        double chi2 = state.getChi2Increment(hit, surface);
        debug() << "    Hit at (" << hit.getPosition()[0] << ", " 
        << hit.getPosition()[1] << ", " << hit.getPosition()[2] 
        << ") has chi2 = " << chi2 << " (limit: " << m_maxChi2 << ")" << endmsg;
        
        // Accept hit if chi-square is below threshold
        if (chi2 < m_maxChi2) {
            compatibleHits.push_back(std::make_tuple(i, hit, surface));
        }
    }
    
    // Sort by chi-square (best hits first)
    std::sort(compatibleHits.begin(), compatibleHits.end(), 
        [&state](const auto& a, const auto& b) {
            return state.getChi2Increment(std::get<1>(a), std::get<2>(a)) < 
                   state.getChi2Increment(std::get<1>(b), std::get<2>(b));
        });
    
    return compatibleHits;
}

TrackState KalmanTracking::predictStep(const TrackState& state, const dd4hep::rec::Surface* nextSurface) {
    return state.predictTo(nextSurface, m_field, *m_materialManager, m_particleProperties);
}

TrackState KalmanTracking::filterStep(const TrackState& predicted, 
                                    const edm4hep::TrackerHitPlane& hit,
                                    const dd4hep::rec::Surface* surface) {
    return predicted.update(hit, surface);
}

void KalmanTracking::createTrack(edm4hep::TrackCollection* trackCollection, const Track& track) const {
    // Create new standard track
    auto std_track = trackCollection->create();
    
    // Set track properties
    std_track.setChi2(track.chi2());
    std_track.setNdf(track.ndf());
    
    // Get final state
    const auto& finalState = track.states().back();
    const auto& params = finalState.parameters();
    
    // Create a new TrackState
    edm4hep::TrackState trackState;
    
    // Convert our (q/pT, phi, eta, d0, z0) parameters to EDM4hep format
    trackState.D0 = params(3);           // d0
    trackState.phi = params(1);           // phi
    trackState.omega = params(0);         // q/pT
    trackState.Z0 = params(4);            // z0
    trackState.tanLambda = std::sinh(params(2)); // Convert eta to tanLambda
    trackState.location = edm4hep::TrackState::AtIP;
    
    // Set reference point
    Eigen::Vector3d refPoint = finalState.positionAtOrigin();
    
    // Use setReferencePoint if available
    float refPointArray[3] = {
        static_cast<float>(refPoint.x()),
        static_cast<float>(refPoint.y()),
        static_cast<float>(refPoint.z())
    };
    trackState.referencePoint = refPointArray;
    
    // Add track state to track
    std_track.addToTrackStates(trackState);
    
    // Add hits to track
    for (const auto& hit : track.hits()) {
        std_track.addToTrackerHits(hit);
    }
}