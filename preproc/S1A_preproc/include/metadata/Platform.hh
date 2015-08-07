#ifndef PLATFORM_HH
#define PLATFORM_HH 1

#include <string>

/**
 * The Platform class holds information about a spaceborne radar satellite
 *
 * @author Walter Szeliga
 */
class Platform {
private:
    std::string mission; //!< The name of the mission
public:
    Platform();
    ~Platform();
    void setMission(const std::string mission);
    std::string getMission() const;
    /*void setPlanet(const Planet planet);
    const Planet getPlanet();*/
};

#endif
