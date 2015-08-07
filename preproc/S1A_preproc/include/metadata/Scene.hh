#ifndef SCENE_HH
#define SCENE_HH 1

#include <string>
#include "boost/date_time/posix_time/posix_time.hpp" 

/**
 * The scene class holds information about a particular data acquisition from a 
 * spaceborne radar satellite 
 *
 * @author Walter Szeliga
 */
class Scene {
private:
    double startingRange; //!< The range from the satellite to the nearest range bin [m]
    double iBias; //!< The in-phase bias of an unfocused pixel
    double qBias; //!< The quadrature bias of an unfocused pixel
    double satelliteHeight; //!< The satellite height measured in the nadir direction to the ground track [m]
    int orbitNumber; //!< The orbit number of the mission
    std::string passDirection; //!< A text string indicating the direction of satellite motion (i.e. Ascending vs. Descending)
    std::string processingFacility; //!< The data processing facility
    std::string processingLevel; //!< The data processing level
    std::string polarization; //!< data polarization
    std::string beam; //!< the beam mode
    int numberOfLines;
    int numberOfSamples;
    boost::posix_time::ptime sensingStart; //!< time in UTC for the first unfocused data line
    boost::posix_time::ptime sensingStop; //!< time in UTC for the last unfocused data line
public:
    Scene();
    ~Scene();
    void setStartingRange(double startingRange);
    void setIBias(double iBias);
    void setQBias(double qBias);
    void setSatelliteHeight(double satelliteHeight);
    void setPassDirection(std::string passDirection);
    void setOrbitNumber(int orbitNumber);
    void setProcessingFacility(std::string processingFacility);
    void setProcessingLevel(std::string processingLevel);
    void setPolarization(std::string polarization);
    void setBeam(std::string beam);
    void setNumberOfLines(int numberOfLines);
    void setNumberOfSamples(int numberOfSamples);
    void setSensingStart(boost::posix_time::ptime sensingStart);
    void setSensingStop(boost::posix_time::ptime sensingStop);
    double getStartingRange() const;
    double getIBias() const;
    double getQBias() const;
    double getSatelliteHeight() const;
    std::string getPassDirection() const;
    int getOrbitNumber() const;
    std::string getProcessingFacility() const;
    std::string getProcessingLevel() const;
    std::string getPolarization() const;
    std::string getBeam() const;
    boost::posix_time::ptime getSensingStart() const;
    boost::posix_time::ptime getSensingStop() const;
    std::string toString();
};

#endif
