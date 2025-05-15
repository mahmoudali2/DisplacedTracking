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
     * Get material parameters at current position
     */
    virtual genfit::Material getMaterialParameters() override {
        // Get material at current position
        dd4hep::Material material = m_manager->materialAt(m_currentPos);
        
        // Convert to GenFit Material
        double density = material.density() / (dd4hep::g/dd4hep::cm3); // Convert to g/cm³
        double Z = material.Z();
        double A = material.A() / (dd4hep::g/dd4hep::mole); // Convert to g/mole
        double radiationLength = material.radLength() / dd4hep::cm; // Convert to cm
        double meanExcitationEnergy = 16.0e-9 * pow(Z, 0.9); // Approximation in GeV
        
        return genfit::Material(density, Z, A, radiationLength, meanExcitationEnergy);
    }
    
    /**
    * Find distance to next material boundary
    * For the moment putting no near boundaries 
    */
    virtual double findNextBoundary(const genfit::RKTrackRep* rep,
                                const genfit::M1x7& state7,
                                double step,
                                bool varField) override {

       // Default large distance (10m)
       return 1.0e6; 
    }
    
private:
    dd4hep::rec::MaterialManager* m_manager;
    dd4hep::rec::Vector3D m_currentPos;
    dd4hep::rec::Vector3D m_currentDir;
};

#endif // GENFIT_ADAPTERS_H