#ifndef SIMULATIONSETTINGS_H
#define SIMULATIONSETTINGS_H

#include "../palabos/src/palabos3D.h"
#include "../palabos/src/palabos3D.hh"   // include full template code

#include <string>
#include <sstream>
#include <chrono>

#include "../palabos/tools/tools.h"
#include "../palabos/tools/tools.hh"
#include "../palabos/tools/plb_physicalFlowParam.h"

using namespace std;
using namespace plb;


template<typename T>
class SimulationSettings {
    public:
        SimulationSettings();
        std::string description;
        std::string outDir;  // Output directory
        std::string histDir, histDirInit; // histogram sub directory
        std::string distDir, distDirInit; // histogram sub directory
        plint numberOfCheckpoints;
        std::string isoSurface, isoSurface_init;
        int STLsmoothing;

        // optimisation
        plint maxSeros;
        char serosConstraint;
        double targetDiameter;
        bool multiLayerSedi;
        bool doSTL;
        bool doHist;
        bool doSedi;
        bool doEros;

        bool stopForKineticEnergy;
        bool stopForDensityAvgDiff;
        bool stopForShearStressDiff;
        T factorSeros;
        T smoothFactor; // exponential smoothing of all measurements
        bool trans; // do temporal averaging of shearStressField between Seros
        T kineticEnergyAvgSmoothGradThreshold;
        T densityAvgDiffSmoothGradSmoothThreshold;
        bool inverseFlow;
        bool openFlow;
        bool freeSlip;

        // geometrical parameters in dimensionless units
        std::string shape;
        char directionOfInlet, directionOfOutlet;
        bool flatMode, freezeInnerSide;
        T xLen,yLen,zLen;
        T inletCentre0,outletCentre0;
        T inletDiameter0,outletDiameter0,inletHeight0,outletHeight0;
        T bufferShiftIn0,bufferShiftOut0;
        T pressureDist0, cross0;
        // geometrical parameters in lattice units
        plint inletCentre;
        plint outletCentre;
        plint inletRadius;
        plint outletRadius;
        plint inletHeight;
        plint outletHeight;
        plint bufferShiftIn;
        plint bufferShiftOut;
        plint cross;
        plint pressureDist;

        plint getOutletCentre(bool left=true);

        Box3D getInletWindow();
        Box3D getInletTwoWindow();
        Box3D getInletRegion();
        Box3D getOutletWindow();
        Box3D getOutletTwoWindow();
        Box3D getMiddleLayer();
        Box3D getBottomLayer();
        Box3D getTopLayer();

        // numerical
        T Re;
        T nodesPerMeter, resolution;
        plint Nx,Ny,Nz;
        T uLB, nuLB, uMax, tau;
        T cSma;
        T uInlet, uDev;
        //IncomprFlowParam<T>* param;
        unique_ptr<PhysicalFlowParam<T>> param;
        T deltaT;
        T deltaX;
        T omega;

        // physical
        T phNu;
        T phL0;
        T phU0;
        T charLength;
        T charVel;
        T charTime;
        T getPhTime(plint it);

        // Iteration intervals for output, checkpointing, etc:
        // variables in physical unit time
        T startTime;
        T checkTime;
        T rampTime;   // time at which upstream is fully developed
        //T rampBeforeSerosTime;   // time at which reshaping may start
        T logTime;
        T analysisTime;
        T delayAnalysisTime;
        T serosTime;
        // variables in iteration steps
        plint charT;
        plint startT;
        plint rampT;
        //plint rampBeforeSerosT;
        plint logT;
        plint analysisT;
        plint delayAnalysisT;
        plint serosT;

        plint imgIntervall;
        plint checkIntervall;
        plint vtkFlowIntervall;
        plint vtkTopoIntervall;
        void read(std::string xmlConfigFile);
        void printWriteLogFile(std::string const& title);
        void createDir();
        string popuName(plint iT);
        string topoName(plint iT);

        // visualisation
        bool tauScaled;
        T rhoMax, tauMax;

        std::string popuNameStem;
        std::string topoNameStem;

        std::string iterFileName;
        std::string serosFileName;
        std::string iterInitializeFileName;
        std::string iterStemFileName;
        std::string serosStemFileName;
};

std::string isoSurfaceCommand(SimulationSettings<T>& sims, bool newBestSeros);


#endif
