#ifndef SIMULATIONSETTINGS_HH
#define SIMULATIONSETTINGS_HH

#include "../palabos/src/palabos3D.h"
#include "../palabos/src/palabos3D.hh"   // include full template code

#include <string>
#include <sstream>
#include <chrono>

#include "../palabos/tools/tools.h"
#include "../palabos/tools/tools.hh"
#include "../palabos/tools/plb_physicalFlowParam.h"
#include "./simulationSettings.h"
#include "./flags.h"

using namespace std;
using namespace plb;


template<typename T>
SimulationSettings<T>::SimulationSettings()
    : iterFileName(  "measurementsFlow0"),
      serosFileName( "measurementsSeros0"),
      iterInitializeFileName(  "measurementsFlow_initialize"),
      iterStemFileName(  "measurementsFlow"),
      serosStemFileName( "measurementsSeros")
{}

template<typename T>
void SimulationSettings<T>::read(std::string xmlConfigFile) {
    XMLreader document(xmlConfigFile);

    document["output"]["description"].read( description);
    document["output"]["img"].read(         imgIntervall );
    document["output"]["vtkFlow"].read(     vtkFlowIntervall );
    document["output"]["vtkTopo"].read(     vtkTopoIntervall );
    document["output"]["check"].read(       checkIntervall );


    document["physical"]["Nu"].read(phNu);
    document["physical"]["L0"].read(phL0);
    document["physical"]["U0"].read(phU0);

    document["optimisation"]["maxSeros"].read(           maxSeros     );
    document["optimisation"]["serosConstraint"].read(    serosConstraint);
    if (serosConstraint == 'D') {
        document["optimisation"]["multiLayerSedi"].read( multiLayerSedi);
        document["optimisation"]["targetDiameter"].read( targetDiameter);
    } else
        multiLayerSedi = false;

    document["optimisation"]["doSTL"].read(doSTL);
    if (doSTL)
        document["optimisation"]["STLsmoothing"].read(STLsmoothing);
    document["optimisation"]["doSedi"].read(doSedi);
    document["optimisation"]["doEros"].read(doEros);
    document["optimisation"]["doHist"].read(doHist);
    document["optimisation"]["factorSeros"].read(factorSeros);

    //http://www.inventoryops.com/articles/exponential_smoothing.htm
    document["optimisation"]["smoothFactor"].read(smoothFactor);
    document["optimisation"]["transient"].read(trans);
    document["optimisation"]["kineticEnergyAvgSmoothGradThreshold"].read(kineticEnergyAvgSmoothGradThreshold);
    document["optimisation"]["densityAvgDiffSmoothGradSmoothThreshold"].read(densityAvgDiffSmoothGradSmoothThreshold);
    document["optimisation"]["stopCriterion"]["kineticEnergy"].read(stopForKineticEnergy);
    document["optimisation"]["stopCriterion"]["densityAvgDiff"].read(stopForDensityAvgDiff);
    document["optimisation"]["stopCriterion"]["shearStressDiff"].read(stopForShearStressDiff);
    document["optimisation"]["inverseFlow"].read( inverseFlow );
    document["optimisation"]["openFlow"].read(    openFlow    );
    document["optimisation"]["freeSlip"].read(    freeSlip    );


    if (Re == 0)
        document["numerics"]["Reynolds"].read(Re);
    if (phL0 == 0) phL0=1.;
    if (phNu == 0) {
        if (phU0 == 0) phU0=1.;
        phNu = phL0*phU0/Re;
    }

    document["numerics"]["uDev"].read(uDev);
    document["numerics"]["nodesPerMeter"].read(nodesPerMeter);
    if (resolution == 0)
        document["numerics"]["resolution"].read(resolution);
    if (resolution == 0)
        resolution = nodesPerMeter * phL0;
    else
        nodesPerMeter = resolution / phL0;

    document["numerics"]["uLB"].read(uLB);
    if (uLB == 0) {
        document["numerics"]["nuLB"].read(nuLB);
        if (nuLB) {
            uLB = nuLB * (T)Re / (T)resolution;
        } else {

            document["numerics"]["uMax"].read(uMax);
            if (uMax == 0) uMax = 0.05; // -> Ma < 0.1
            document["numerics"]["tau"].read(tau);
            if (tau == 0) tau = 0.8;    // recommended for low Re
            nuLB = (tau - 0.5) / (T)3;  // for tau=0.8 -> nu=0.1
            uLB = nuLB * (T)Re / (T)resolution;
            if (uLB > uMax) // happens at high Re
                // reduce nu (or increase resolution)
                uLB = uMax;
                // will lead to nuLB = uMax * (T)resolution / (T)Re;
        }
    }

    document["numerics"]["cSma"].read(cSma);

    PLB_ASSERT(util::greaterThan_abs(uLB, (T) 0));

    document["time"]["start"].read(   startTime   );
    document["time"]["ramp"].read(    rampTime    );
    document["time"]["log"].read(     logTime     );
    document["time"]["analysis"].read( analysisTime );
    document["time"]["delayAnalysis"].read( delayAnalysisTime );
    document["time"]["seros"].read(    serosTime    );


    document["geometry"]["shape"].read( shape );
    if (shape != "L" && shape != "T" && shape != "U" && shape != "-" &&
        shape != "~" && shape != "D" && shape != "Y" && shape != "Ti" && shape != "d") {
        pcout << "Only \"L\", \"T\", \"Ti\", \"U\", \"-\", \"~\", \"D\" & \"Y\" -shape supported for now! Sorry." << endl;
        exit(1);
    }
    if (shape == "~" || shape == "-" || shape == "D") uInlet = uLB;
    else uInlet = -1*uLB;

    if (shape == "~" || shape == "-" || shape == "D")
        directionOfInlet = 'x';
    else
        directionOfInlet = 'y';
    if (shape == "L" || shape == "Ti" ||shape == "T" )
        directionOfOutlet = 'x';
    else
        directionOfOutlet = 'y';


    document["geometry"]["flatMode"].read( flatMode );
    document["geometry"]["xLen"].read( xLen );
    document["geometry"]["yLen"].read( yLen );
    document["geometry"]["zLen"].read( zLen );
    // the following will result in Nz = minLayer+1
    T minLayer = 3.;
    if( zLen == -1 ) zLen = minLayer/resolution;
    document["geometry"]["inletPosition"].read(   inletCentre0   );
    document["geometry"]["outletPosition"].read(  outletCentre0  );
    document["geometry"]["inletDiameter"].read(   inletDiameter0   );
    document["geometry"]["outletDiameter"].read(  outletDiameter0  );
    document["geometry"]["inletHeight"].read( inletHeight0 );
    document["geometry"]["outletHeight"].read( outletHeight0 );
    document["geometry"]["bufferShiftIn"].read(  bufferShiftIn0     );
    document["geometry"]["bufferShiftOut"].read( bufferShiftOut0    );
    document["geometry"]["cross"].read( cross0 );
    document["geometry"]["pressureDist"].read( pressureDist0 );
    document["geometry"]["freezeInnerSide"].read( freezeInnerSide );

    inletCentre     = resolution*inletCentre0;
    outletCentre    = resolution*outletCentre0;
    inletRadius     = resolution*inletDiameter0/2;
    outletRadius    = resolution*outletDiameter0/2;
    inletHeight    = resolution*inletHeight0;
    outletHeight    = resolution*outletHeight0;
    bufferShiftIn   = resolution*bufferShiftIn0;
    bufferShiftOut  = resolution*bufferShiftOut0;
    cross           = resolution*cross0;
    pressureDist    = resolution*pressureDist0;


    document["visualisation"]["rhoMax"].read(rhoMax);
    document["visualisation"]["tauScaled"].read(tauScaled);
    document["visualisation"]["tauMax"].read(tauMax);

    /*
    Numeric = new IncomprFlowParam<T>(
            (T) 1.,     // physicalU
            (T) uLB,   // latticeU
            (T) Re,
            (T) outletDiameter,     // physicalLength
            resolution,
            xLen,
            yLen,
            zLen
    );
    */

    param = make_unique<PhysicalFlowParam<T>>(PhysicalFlowParam<T>(
            uLB,            // latticeU
            Re,              // Reynolds number
            phNu,            // physical kinematic viscosity
            phL0,            // physicalLength
            nodesPerMeter,
            xLen,
            yLen,
            zLen
    ));


    Nx = param->getNx();
    Ny = param->getNy();
    Nz = param->getNz();
    /*
    if (Nx%2==0 || Ny%2==0) {
        pcout << "Nx and Ny are not odd numbers, this might result in undefined behavior!" << endl;
        pcout << "Please use other nodesPerMeter/resolution and try again." << endl;
        exit(1);
    }
    */
    deltaT = param->getDeltaT();
    deltaX = param->getDeltaX();
    omega = param->getOmega();
    charTime   = 1.;

    charT         = param->nStep(charTime);
    startT        = param->nStep(startTime);
    rampT         = param->nStep(rampTime);
    //rampBeforeSerosT = param->nStep(rampBeforeSerosTime);
    logT          = param->nStep(logTime);
    analysisT      = param->nStep(analysisTime);
    delayAnalysisT = param->nStep(delayAnalysisT);
    serosT        = param->nStep(serosTime);

    int dimensions;
    if (flatMode) dimensions = 2;
    else          dimensions = 3;

    ostringstream sOut;  // 2D_L_Re_x_description
    sOut << dimensions << "D_" << shape << "_Re" << Re << "_x" << resolution
         << "_s" << serosTime << "_r" << factorSeros;
    if (serosConstraint == 'D' || serosConstraint == 'V' ) {
        sOut << "_" << serosConstraint;
        if (serosConstraint == 'D' && !multiLayerSedi)
            sOut << "s";
    }
    if (trans)
        sOut << "_Trans";
    if (inverseFlow)
        sOut << "_inverse";
    if (openFlow)
        sOut << "_open";
    if (description != "")
        sOut << "_" << description;
    outDir = sOut.str()+'/';;
    global::directories().setOutputDir(outDir);

    histDir = outDir+"histogram";
    histDirInit = histDir+"-init";
    distDir = outDir+"distribution";
    distDirInit = distDir+"-init";
    isoSurface = outDir+"isoSurface.ply";
    isoSurface_init = outDir+"isoSurface_init.ply";

    ostringstream sPopu; // popu2D_L_Re_x_Nx_Ny_Nz
    sPopu << "popu" << dimensions << "D_" << shape << "_Re" << Re << "_x" << resolution << "_Nx" << Nx << "_Ny" << Ny << "_Nz" << Nz;
    popuNameStem = sPopu.str();

    ostringstream sTopo; // topo2D_L_x_Nx_Ny_Nz
    sTopo << "topo" << dimensions << "D_" << shape <<                "_x" << resolution << "_Nx" << Nx << "_Ny" << Ny << "_Nz" << Nz;
    topoNameStem = sTopo.str();

}

template<typename T>
plint SimulationSettings<T>::getOutletCentre(bool left) {
    if ( shape == "Y" ) {
        if (left)
            return Nx/2 - outletCentre;
        else
            return Nx/2 + outletCentre;
    } else {
        return outletCentre; }
}

template<typename T>
Box3D SimulationSettings<T>::getInletWindow() {
    if ( shape == "~" || shape == "-" || shape == "D")
        return Box3D(0,0, inletCentre-inletRadius,inletCentre+inletRadius, 0,Nz-1);
    else
        return Box3D(inletCentre-inletRadius,inletCentre+inletRadius, Ny-1,Ny-1, 0,Nz-1); }
template<typename T>
Box3D SimulationSettings<T>::getInletTwoWindow() {
    return Box3D(inletCentre-inletRadius,inletCentre+inletRadius, 0,0, 0,Nz-1); } // only to be used with Ti-shape
template<typename T>
Box3D SimulationSettings<T>::getInletRegion() {
    if (shape == "~" || shape == "-" || shape == "D")
        return Box3D(pressureDist,pressureDist+2, inletCentre-inletRadius,inletCentre+inletRadius, 0,Nz-1);
    else
        return Box3D(inletCentre-inletRadius,inletCentre+inletRadius, Ny-1-pressureDist-2,Ny-1-pressureDist, 0,Nz-1); }
template<typename T>
Box3D SimulationSettings<T>::getOutletWindow() {
    switch (shape.c_str()[0]) {
        case 'U': return Box3D(outletCentre-outletRadius,outletCentre+outletRadius, Ny-1,Ny-1, 0,Nz-1);
        case 'Y': return Box3D(Nx/2-outletCentre-outletRadius,Nx/2-outletCentre+outletRadius, 0,0, 0,Nz-1);
        case '-': return Box3D(Nx-1,Nx-1, inletCentre-inletRadius,inletCentre+inletRadius, 0,Nz-1);
                  // because for '-' shape, outlet is on same height as inlet;
                  // outletCentre is used for position of diversion
        default : return Box3D(Nx-1,Nx-1, outletCentre-outletRadius,outletCentre+outletRadius, 0,Nz-1); } }
template<typename T>
Box3D SimulationSettings<T>::getOutletTwoWindow() {
    switch (shape.c_str()[0]) {
        case 'Y': return Box3D(Nx/2+outletCentre-outletRadius,Nx/2+outletCentre+outletRadius, 0,0, 0,Nz-1);
        default : return Box3D(0,0, outletCentre-outletRadius,outletCentre+outletRadius, 0,Nz-1); } } // only to be used with T-shape
template<typename T>
Box3D SimulationSettings<T>::getMiddleLayer() {
    return Box3D(0, Nx-1, 0, Ny-1, Nz/2, Nz/2); }
template<typename T>
Box3D SimulationSettings<T>::getBottomLayer() {
    return Box3D(0, Nx-1, 0, Ny-1, 0, 0); }
template<typename T>
Box3D SimulationSettings<T>::getTopLayer() {
    return Box3D(0, Nx-1, 0, Ny-1, Nz-1, Nz-1); }

template<typename T>
T SimulationSettings<T>::getPhTime(plint iT) {
    return param->getPhysicalT(iT*deltaT); }

template<typename T>
void SimulationSettings<T>::printWriteLogFile(std::string const& title) {
    time_t timer;
    char timebuffer[26];
    struct tm* tm_info;
    time(&timer);
    tm_info = localtime(&timer);
    strftime(timebuffer, 26, "%Y-%m-%d_%H:%M:%S", tm_info);

    stringstream logText;
    logText << timebuffer << endl
    << title << endl << endl
    << "--------------- Simulation Parameters:-----------" << endl
    << " Start   (t) = " << startTime   << endl
    << "        (iT) = " << startT      << endl
    << " End (seros) = " << maxSeros    << endl
    << " InletDiameter:  Ni = " << inletRadius*2+1         << endl
    << " OutletDiameter: No = " << outletRadius*2+1        << endl;
    if (inletHeight != 0)
        logText << " InletHeight: NoH = " << inletHeight+1 << endl;
    if (outletHeight != 0)
        logText << " OutletHeight: NoH = " << outletHeight+1 << endl;
    logText
    << "--------------- Seros Parameters: ---------------" << endl
    << " Velocity (lattice units): uLB = " << uLB  << endl;
    if (uDev != 0) {
        logText << "              deviation:  uDev = " << uDev << endl;
    }
    logText
    << " Initialization time:      (t) = " << rampTime     << endl
    << "                          (iT) = " << rampT        << endl
    << " Seros intervall:          (t) = " << serosTime    << endl
    << "                          (iT) = " << serosT       << endl
    << " Flow analysis intervall:  (t) = " << analysisTime  << endl
    << "                          (iT) = " << analysisT     << endl
    << " delay analysis:           (t) = " << delayAnalysisTime << endl
    << "                          (iT) = " << delayAnalysisT    << endl
    << " temporal averaging            = " << trans    << endl
    << " Seros constraint: constant " << serosConstraint   << endl;
    if (serosConstraint == 'D')
        logText << " add multiple layers during sedi (1=yes/0=no): " << multiLayerSedi << endl;
    logText
    << "-------------------------------------------------" << endl;
    // print to screen
    pcout << logText.str();
    // save as file
    rewriteFile(logText.str(), outDir+"serosLog.dat");

    // do the same for PhysicalFlowParam
    writeLogFile(*param, outDir+"plbLog.dat");
}


template<typename T>
string SimulationSettings<T>::popuName(plint iT) {
    if (iT == 0) return popuNameStem;
    stringstream sPopuName;
    sPopuName << popuNameStem
              << "_iT" << iT << "_time" << param->stepToDimless(iT);
    return sPopuName.str();
}

template<typename T>
string SimulationSettings<T>::topoName(plint iT) {
    if (iT == 0) return topoNameStem;

    stringstream sTopoName;
    sTopoName << topoNameStem;
    if (iT == -1)
        sTopoName << "_latestSeros";
    else if (iT == -2)
        sTopoName << "_bestSeros";
    else
        sTopoName << "_iT" << iT << "_time" << param->stepToDimless(iT);
    return sTopoName.str();
}

std::string isoSurfaceCommand(SimulationSettings<T>& sims, bool newBestSeros)
{
    int contour = flagBB;
    stringstream command;
    std::string executable = "./../palabos/externalLibraries/IsoSurface/Bin/Linux/IsoSurfaceExtraction";
    std::string topoBinaryLocation = sims.outDir+sims.topoName(-1)+".binary";
    std::string extractorLog = sims.outDir+"extractorLog";
    command << executable << " --in " << topoBinaryLocation
            << " --res " << sims.Nx <<" " << sims.Ny <<" " << sims.Nz
            << " --iso " << contour << " --sIters " << sims.STLsmoothing;
    if (newBestSeros)
        command << " --out " << sims.outDir << "isoSurface_best.ply";
    else
        command << " --out " << sims.isoSurface;
    command <<" >> " << extractorLog;

    return command.str();
}


#endif
