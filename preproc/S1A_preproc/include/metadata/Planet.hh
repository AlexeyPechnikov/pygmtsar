#ifndef PLANET_HH
#define PLANET_HH 1

/**
 * The Planet class holds information about planetary constants
 *
 * @author Walter Szeliga
 */
class Planet {
public:
    static const Planet Earth; //< The Earth
private:
    double mass; //!< The mass of the planet [kg]
    double spinRate; //!< The rotation rate [rad/s]
    double semiMajorAxis; //!< The radius of the semi-major axis [m]
    double semiMinorAxis; //!< The radius of the semi-minor axis [m]
    double eccentricity; //!< The eccentricity of the radius
    double flattening; //!< The flattening

    Planet(double mass, double spinRate, double semiMajorAxis, double semiMinorAxis, double eccentricity, double flattening) {
        this->mass = mass;
        this->spinRate = spinRate;
        this->semiMajorAxis = semiMajorAxis;
        this->semiMinorAxis = semiMinorAxis;
        this->eccentricity = eccentricity;
        this->flattening = flattening;
    }
public:
    const double getMass() const;
    const double getSpinRate() const;
    const double getSemiMajorAxis() const;
    const double getSemiMinorAxis() const;
    const double getEccentricity() const;
    const double getFlattening() const;
    const double getGM() const;
    const double getEccentricitySquared() const;
};

#endif
