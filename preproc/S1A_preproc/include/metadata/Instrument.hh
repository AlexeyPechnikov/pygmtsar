#ifndef INSTRUMENT_HH
#define INSTRUMENT_HH 1

#include <string>
#include "Platform.hh"

/**
 * The Instrument class holds information about a radar instrument on-board
 * a particular spaceborne satellite platform
 *
 * @author Walter Szeliga
 */
class Instrument {
private:
    Platform *platform; //!< The platform the instrument is attached to
    double wavelength; //!< The radar wavelength [m]
    double incidenceAngle; //!< The incidence angle to the firt range bin
    double pulseRepetitionFrequency; //!< The pulse repetition frequency [Hz]
    double rangePixelSize; //<! The range pixel size [m]
    double azimuthPixelSize; //<! The azimuth pixel size [m]
    double pulseLength; //!< The pulse length [s]
    double chirpSlope; //!< The chirp slope [Hz/s]
    double antennaLength; //!< The length of the antenna [m]
    int antennaSide; //!< An integer representing the look direction of the satellite
    double rangeSamplingFrequency; //!< The range sampling frequency [Hz]
public:
    Instrument();
    ~Instrument();
    void setPlatform(Platform *platform); 
    void setWavelength(double wavelength); 
    void setIncidenceAngle(double incidenceAngle);
    void setPulseRepetitionFrequency(double pulseRepetitionFrequency);
    void setRangePixelSize(double rangePixelSize);
    void setAzimuthPixelSize(double azimuthPixelSize);
    void setPulseLength(double pulseLength);
    void setChirpSlope(double chirpSlope);
    void setAntennaLength(double antennaLength);
    void setAntennaSide(int antennaSide);
    void setRangeSamplingFrequency(double rangeSamplingFrequency);
    Platform *getPlatform();
    double getWavelength() const;
    double getIncidenceAngle() const;
    double getPulseRepetitionFrequency() const;
    double getRangePixelSize() const;
    double getAzimuthPixelSize() const;
    double getPulseLength() const;
    double getChirpSlope() const;
    double getAntennaLength() const;
    int getAntennaSide() const;
    double getRangeSamplingFrequency() const;
    std::string toString();
};

#endif
