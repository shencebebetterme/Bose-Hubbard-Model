#include "pch.h"

extern int numSites;
extern int numParticles;
extern double intStrength; //the interaction strength U/2 (with J set to 1)

static enum argStr
{
	num_Sites,
	num_Particles,
	inter_Strength
};

using mapType = std::map<std::string, argStr>;

void showHelp();
void initMap(mapType&);
void parseAllArgs(mapType&, int, char**);

void argParser(int argc, char* argv[]) {
	//no arguments, use the default value
	if (argc == 1) {
		return;
	}

	//show help message
	if (std::string(argv[1]) == "-h") {
		showHelp();
	}

	mapType argMap;
	initMap(argMap);
	parseAllArgs(argMap, argc, argv);
}

// initialize the map
void initMap(mapType& map) {
	map["-nS"] = num_Sites;
	map["-nP"] = num_Particles;
}

// show help message
void inline showHelp() {
	std::cout << "\npass command line arguments\n"
		<< "-nS\tnumber of Sites, default value = 1\n"
		<< "-nP\tnumber of Particles, default value = 1\n"
		<< "-i\trelative interaction strength U/2J, default value = 1\n"
		<< std::endl;
	std::exit(EXIT_FAILURE);
}

// pass all arguments (starting from the second arg)
void inline parseAllArgs(mapType& map, int argc, char* argv[]) {
	for (int i = 1; i < argc; i++) {
		std::string arg_i(argv[i]);
		// other command arguments have structure arg=val
		int equalSign = arg_i.find("=");
		// obtain the arg-val pair
		std::string arg_ic = arg_i.substr(0, equalSign);//arg_i content
		std::string arg_iv = arg_i.substr(equalSign + 1, -1);//arg_i value
		switch (map[arg_ic])
		{
		case num_Sites: numSites = std::stoi(arg_iv); break;
		case num_Particles: numParticles = std::stoi(arg_iv); break;
		}
	}
}