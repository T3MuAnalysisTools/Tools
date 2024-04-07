#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

// Define a functor class
class Track_Functor : public ROOT::Math::IBaseFunctionMultiDim {
public:
    // Constructor with fixed constants
    Track_Functor(TrackParticle particle1, TrackParticle particle2) : m_particle1(particle1), m_particle2(particle2) {}
    
    // Define the distance function
    double DoEval(const double *params) const override {
        double t1 = params[0];
        double t2 = params[1];

        TVector3 helix1 = computeHelixCoordinates(m_particle1, t1);
        TVector3 helix2 = computeHelixCoordinates(m_particle2, t2);

        // Compute the Euclidean distance between the points
        double distance = (helix1-helix2).Mag();

        return distance;
    }
    
    // Clone method for deep copying
    virtual Track_Functor* Clone() const override {
        return new Track_Functor(*this);
    }

    // Function dimensionality
    unsigned int NDim() const override {
        return 2; // The function is two-dimensional
    }

private:
    TrackParticle m_particle1;
    TrackParticle m_particle2;
    
    TVector3 computeHelixCoordinates(const TrackParticle& particle, double &z) const {
        double kappa=particle.Parameter(TrackParticle::kappa);
        double lambda=particle.Parameter(TrackParticle::lambda);
        double phi0=particle.Parameter(TrackParticle::phi);
        double dxy=particle.Parameter(TrackParticle::dxy);
        double dz=particle.Parameter(TrackParticle::dz);

        double s=(z-dz)/tan(lambda);
        //double radius=2.0/kappa; // Radius of curvature
        //double radius=1.0/(kappa*particle.BField()); // Radius of curvature
        double radius=1.0/(kappa); // Radius of curvature
        
        // Calculate the helix parameters
        double x=radius*sin(2.0*s*kappa+phi0)-(radius+dxy)*sin(phi0);
        double y=-radius*cos(2.0*s*kappa+phi0)+(radius+dxy)*cos(phi0);
        
        //std::cout << "del_phi: " << del_phi <<", x: " << x << ", y " << y << " and z: " << z << std::endl;
        
        return TVector3(x,y,z);
    }
    
};

