#ifndef GENFIT_ADAPTERS_H
#define GENFIT_ADAPTERS_H

#include "AbsBField.h"
#include "AbsMaterialInterface.h"
#include "DD4hep/Fields.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Vector3D.h"
#include "TVector3.h"

/**
 * Adapter class to make DD4hep magnetic field usable in GenFit
 */
class DD4hepFieldAdapter : public genfit::AbsBField {
public:
    DD4hepFieldAdapter(const dd4hep::OverlayedField& field) : m_field(field) {}
    
    /**
     * Get magnetic field at position (in cm)
     */
    TVector3 get(const TVector3& pos) const override {
        // Convert TVector3 to dd4hep::Position (pos is in cm as required by DD4hep)
        dd4hep::Position ddPos(pos.X(), pos.Y(), pos.Z());
        
        // Get field at position (DD4hep returns field in tesla)
        dd4hep::Direction fieldVector = m_field.magneticField(ddPos);
        
        // Convert to TVector3 (units are kGauss for GenFit, 1 Tesla = 10 kGauss)
        return TVector3(fieldVector.x() * 10.0, 
                        fieldVector.y() * 10.0, 
                        fieldVector.z() * 10.0);
    }
    
private:
    dd4hep::OverlayedField m_field;
};

/**
 * Adapter class to make DD4hep material information usable in GenFit
 */
class DD4hepMaterialAdapter : public genfit::AbsMaterialInterface {
public:
    DD4hepMaterialAdapter(dd4hep::rec::MaterialManager* manager) : m_manager(manager) {}
    
    /**
     * Initialize track parameters - matches the exact signature from AbsMaterialInterface
     */
    virtual bool initTrack(double posX, double posY, double posZ,
                          double dirX, double dirY, double dirZ) override {
        // Store the initial position for later use
        m_currentPos = dd4hep::rec::Vector3D(posX, posY, posZ);
        m_currentDir = dd4hep::rec::Vector3D(dirX, dirY, dirZ);
        
        return true;
    }
    
    /**
     * Get material parameters at current position.
     * Returns vacuum-equivalent material when the position is outside the
     * world volume (materialAt throws in that case).
     */
    virtual genfit::Material getMaterialParameters() override {
        try {
            dd4hep::Material material = m_manager->materialAt(m_currentPos);
            double density  = material.density() / (dd4hep::g/dd4hep::cm3);
            double Z        = material.Z();
            double A        = material.A() / (dd4hep::g/dd4hep::mole);
            double radLen   = material.radLength() / dd4hep::cm;
            double meanExcE = 16.0e-9 * std::pow(Z, 0.9);  // approximate, in GeV
            return genfit::Material(density, Z, A, radLen, meanExcE);
        } catch (...) {
            // Position is outside the geometry world — return vacuum
            return genfit::Material(0.0, 0.0, 1.0, 1e9, 1e-9);
        }
    }
    
    /**
     * Find distance to next material boundary.
     *
     * Returns a step limit that prevents GenFit's RK integrator from jumping
     * over a material boundary.  We use a conservative fixed step of 0.5 cm
     * in the muon system — small enough to capture boundaries but large enough
     * to keep the step count reasonable.  A full boundary-search would need
     * per-step geometry traversal that is expensive; for the muon system the
     * simple cap is adequate.
     */
    virtual double findNextBoundary(const genfit::RKTrackRep* /*rep*/,
                                const genfit::M1x7& state7,
                                double step,
                                bool /*varField*/) override {
        // Clamp to a conservative maximum step
        // state7 = {x, y, z, px/p, py/p, pz/p, q/p} in cm
        const double maxStep = 0.5;  // cm — safe for the muon detector geometry
        return std::min(std::abs(step), maxStep);
    }
    
private:
    dd4hep::rec::MaterialManager* m_manager;
    dd4hep::rec::Vector3D m_currentPos;
    dd4hep::rec::Vector3D m_currentDir;
};

#endif // GENFIT_ADAPTERS_H