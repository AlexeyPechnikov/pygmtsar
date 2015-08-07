/*
 *  Copyright (C) 2012 Walter M. Szeliga
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <string>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include "csk.hh"
#include "utilities.hh"
#include "gmtsar.hh"
#include "Doppler.hh"
#include "config.h"

int
main(int argc, char *argv[])
{
	std::string filename;
	std::string prefix;
	namespace po = boost::program_options;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h","print help message")
		("version,v'","print version number")
		("prefix,p",po::value<std::string>(),"output file prefix")
		("input,i",po::value<std::string>(),"input file name");

	po::variables_map vm;
	try {
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);
	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}

	if (vm.count("version")) {
		std::cout << "make_raw_csk " << PACKAGE_VERSION << "\n"
			<< "Copyright (C) 2012 Walter Szeliga\n"
			<< "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.htm>\n"
			<< "This is free software: you are free to change and redistribute it.\n"
			<< "There is NO WARRANTY, to the extent permitted by law" << std::endl;
		return 0;
	}

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		std::cout << "Report bugs to: walter <at> geology.cwu.edu" << std::endl;
    		return 1;
	}

	if (vm.count("input"))
	{
		filename = vm["input"].as<std::string>();
	}
	else
	{
		std::cerr << "input file required" << std::endl;
		return 1;
	}

	if (vm.count("prefix"))
	{
		prefix = vm["prefix"].as<std::string>();
	}
	else
	{
		std::cerr << "output file prefix required" << std::endl;
		return 1;
	}

	// Cook up the output file names
	std::string outputImage = prefix + ".raw";
        std::string outputOrbit = prefix + ".LED";
        std::string outputRSC = prefix + ".PRM";

	CSK csk = CSK(filename);
	std::cout << "Extracting Image" << std::endl;
	SingleBandImage<std::complex<unsigned char> > *image = csk.extractRawImage(outputImage);
	if (image == NULL)
	{
		return 1;
	}

	GMTSAR::writeEphemeris(outputOrbit,csk.getOrbit(),csk.getScene());
	GMTSAR::writeResourceFile(outputRSC,csk.getScene(),csk.getInstrument(),csk.getOrbit(),image);

	delete image;
	return 0;
}
