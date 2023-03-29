#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include "../palabos/src/palabos3D.h"
#include "../palabos/src/palabos3D.hh"   // include full template code
#include "../palabos/externalLibraries/easyppm/easyppm.h"
#include "../palabos/externalLibraries/easyppm/easyppm.c"

#include <string>
#include <vector>

#include "../palabos/tools/tools_geometric.h"
#include "../palabos/tools/tools_geometric.hh"
#include "../palabos/tools/tools_plb.h"
#include "../palabos/tools/tools_plb.hh"
#include "./flags.h"
#include "./geometry.h"

using namespace plb;
using namespace plb::descriptors;
using namespace std;


template<typename T>
Geometry<T>::Geometry(MultiBlockLattice3D<T, NSDESCRIPTOR>& lattice_,
        bool flatMode_)
    : lattice(lattice_),
      Nx(lattice.getNx()),
      Ny(lattice.getNy()),
      Nz(lattice.getNz()),
      topology(
          MultiScalarField3D<int>(
              MultiBlockManagement3D (
                  SparseBlockStructure3D(
                      createRegularDistribution3D(Nx,Ny,Nz)),
                  defaultMultiBlockPolicy3D().getThreadAttribution(),
                  2),
              defaultMultiBlockPolicy3D().getBlockCommunicator(),
              defaultMultiBlockPolicy3D().getCombinedStatistics(),
              defaultMultiBlockPolicy3D().getMultiScalarAccess<int>())
      ),
      topologyForSTL( MultiScalarField3D<int>(Nx,Ny,Nz,0)),
      eroding(        MultiScalarField3D<int>(Nx,Ny,Nz,0)),
      flatMode(flatMode_)
{ }

template<typename T>
void Geometry<T>::updateTopology(MultiScalarField3D<int>& topo, bool updateDynamics) {
    MaskedReplaceAlphaFunctional3D<int, int> removeFlagEroding(flagWet, flagEroding);
        applyProcessingFunctional(
            removeFlagEroding, lattice.getBoundingBox(), topo, topo);
    if (flatMode) {
        // only calculate update in bottom layer:
        applyProcessingFunctional(
            new UpdateTopologyFunctional3D<int>, layer(Nx,Ny,0), topo);
        // and then copy topology to all other layers:
        for(plint iZ=0; iZ<Nz; iZ++)
            copy(topo, layer(Nx,Ny,0), topo, layer(Nx,Ny,iZ));
    } else {
        applyProcessingFunctional(
            new UpdateTopologyFunctional3D<int>, lattice.getBoundingBox(), topo);
    }
    if (updateDynamics) {
        defineDynamics(lattice, topo,
            new BounceBack<T,NSDESCRIPTOR>(0.), flagSolid);
        defineDynamics(lattice, topo,
            new BounceBack<T,NSDESCRIPTOR>(1.), flagSedimenting);
    }
    MaskedReplaceAlphaFunctional3D<int, int> removeFlagSedimenting(flagBB, flagSedimenting);
    applyProcessingFunctional(
        removeFlagSedimenting, lattice.getBoundingBox(), topo, topo);
}

template<typename T>
void Geometry<T>::createTopology(std::string shape, plint inletCentre, plint inletRadius,
        plint outletCentre, plint outletRadius, plint bufferShiftIn, plint bufferShiftOut,
        plint cross, bool freeze, bool openFlow, plint inletHeight, plint outletHeight) {

    if (flatMode) {
        // only calculate topology in bottom layer:
        applyProcessingFunctional(
            new CreateTopologyFunctional3D(
                shape, Nx, Ny, 0,
                inletCentre, inletRadius, outletCentre, outletRadius,
                bufferShiftIn, bufferShiftOut, cross,
                openFlow, inletHeight, outletHeight),
            layer(Nx,Ny,0), topology);
        // and then copy topology to all other layers:
        for(plint iZ=0; iZ<Nz; iZ++)
            copy(topology, layer(Nx,Ny,0), topology, layer(Nx,Ny,iZ));
    } else {
        applyProcessingFunctional(
            new CreateTopologyFunctional3D(
                shape, Nx, Ny, Nz,
                inletCentre, inletRadius, outletCentre, outletRadius,
                bufferShiftIn, bufferShiftOut, cross,
                openFlow, inletHeight, outletHeight),
            lattice.getBoundingBox(), topology);
    }

    if (shape=="U" && freeze) {
        Box3D freezeThisVolume(inletCentre+inletRadius*.5,outletCentre-outletRadius*.5, (Ny-1)/4+outletRadius*.5,Ny-1, 0,Nz-1);
        applyProcessingFunctional(
            new FreezeTopologyFunctional3D<int>(), freezeThisVolume, topology);
    }
    updateTopology(topology);
}

template<typename T>
void Geometry<T>::loadTopologyFromPPM(std::string ppmName) {
    if (flatMode)
        applyProcessingFunctional(
            new PPM2TopologyFunctional3D<int>(ppmName),
            lattice.getBoundingBox(), topology);
}

template<typename T>
void Geometry<T>::loadTopology(std::string fileName) {
    plb_ifstream ifile(fileName.c_str());
    ifile >> topology;
}

template<typename T>
void Geometry<T>::saveTopology(std::string fileName, bool binary) {
    if (binary) {
        saveBinaryBlock(topology, fileName+".binary");
    } else {
        plb_ofstream ofile((fileName).c_str());
        ofile << topology;
    }
}

template<typename T>
void Geometry<T>::markErodingCellsInTopology(MultiScalarField3D<int> &topo) {
    MaskedReplaceAlphaFunctional3D<int, int> markEroding(flagEroding, 1);
    applyProcessingFunctional(
        markEroding, lattice.getBoundingBox(), topo, eroding);
}

template<typename T>
void Geometry<T>::erode(MultiScalarField3D<T> &shearStressNormField, T erosThreshold) {

    // only look at nodes which are marked as flagWet
    // if the shearStress in those nodes is above threshold -> mark as 'eroding'
    if(flatMode) {
        // only calculate erosion in bottom layer:
        applyProcessingFunctional(
            new ErodeFunctional3D<int,T>(erosThreshold),
                layer(Nx,Ny,0), topology, shearStressNormField);
        // and then copy topology to all other layers:
        for(plint iZ=1; iZ<Nz; iZ++)
            copy(topology, layer(Nx,Ny,0), topology, layer(Nx,Ny,iZ));
    } else {
        applyProcessingFunctional(
            new ErodeFunctional3D<int,T>(erosThreshold),
                lattice.getBoundingBox(), topology, shearStressNormField);
    }

    // but they still have to actually erode (set Dynamics before next c&s!)

    // track cells which became fluid cells due to erosion,
    // so they cannot sedi during ongoing reshaping
    eroding = MultiScalarField3D<int>(Nx,Ny,Nz,0);
    MaskedReplaceAlphaFunctional3D<int, int> trackEroding(1, flagEroding);
    applyProcessingFunctional(
        trackEroding, lattice.getBoundingBox(), eroding, topology);

    // recently eroded nodes should not be able to sediment in next step
    // will be marked as flagWet in next UpdateTopologyFunctional
}

template<typename T>
void Geometry<T>::sediment(MultiScalarField3D<int> &topo,
                           MultiScalarField3D<T> &shearStressNormField,
                           T sediThreshold) {

    // only look at nodes which are marked as flagWet (did not erode/sediment in last step)
    // if the shearStress in those nodes is below threshold -> mark as 'sedimenting'
    if (flatMode) {
        // only calculate sedimentation in bottom layer:
        applyProcessingFunctional(
            new SedimentFunctional3D<int,T>(sediThreshold), layer(Nx,Ny,0), topo, shearStressNormField);
        // and then copy topology to all other layers:
        for(plint iZ=0; iZ<Nz; iZ++)
            copy(topo, layer(Nx,Ny,0), topo, layer(Nx,Ny,iZ));
    } else {
        applyProcessingFunctional(
            new SedimentFunctional3D<int,T>(sediThreshold), lattice.getBoundingBox(), topo, shearStressNormField);
    }

    // but they still have to sediment. this will be done in "updateTopology"
}

template<typename T>
MultiScalarField3D<int>& Geometry<T>::getTopology() {
    return topology;
}

template<typename T>
MultiScalarField3D<int>& Geometry<T>::getTopologyForSTL() {
    topologyForSTL = MultiScalarField3D<int>(topology);
    MaskedReplaceAlphaFunctional3D<int, int> markConstraintBB(flagBB, flagConstraintBB);
    MaskedReplaceAlphaFunctional3D<int, int> markBulk(flagWet, flagBulk);
    MaskedReplaceAlphaFunctional3D<int, int> markBuffer(flagWet, flagBuffer);
    applyProcessingFunctional(
        markConstraintBB, lattice.getBoundingBox(), topologyForSTL, topology);
    applyProcessingFunctional(
        markBulk, lattice.getBoundingBox(), topologyForSTL, topology);
    applyProcessingFunctional(
        markBuffer, lattice.getBoundingBox(), topologyForSTL, topology);
    return topologyForSTL;

}


/* ******** CreateTopologyFunctional3D ************************************** */
CreateTopologyFunctional3D::CreateTopologyFunctional3D(
        std::string shape_, plint Nx, plint Ny, plint Nz,
        plint inletCentre,  plint inletRadius, plint outletCentre, plint outletRadius,
        plint bufferShiftIn_, plint bufferShiftOut_, plint cross_,
        bool openFlow_, plint inletHeight, plint outletHeight)
    : shape(shape_),
      X1(Nx-1),
      Y1(Ny-1),
      Z1(Nz-1),
      InC    ( inletCentre               ),
      InR    ( inletRadius               ),
      InRSq  ( inletRadius*inletRadius   ),
      InR2   ( inletHeight/2             ),
      InR2Sq ( InR2*InR2                 ),
      OutC   ( outletCentre              ),
      OutR   ( outletRadius              ),
      OutRSq ( outletRadius*outletRadius ),
      OutR2  ( outletHeight/2            ),
      OutR2Sq( OutR2*OutR2               ),
      bufferShiftIn   ( bufferShiftIn_ ),
      bufferShiftOut  ( bufferShiftOut_ ),
      cross     ( cross_    ),
      openFlow  ( openFlow_ )
{ }


void CreateTopologyFunctional3D::process (
        Box3D domain, ScalarField3D<int>& topology )
{
    // Access the position of the scalarfield
    Dot3D relativePosition = topology.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Convert local coordinates to global ones.
                plint gX = iX + relativePosition.x;
                plint gY = iY + relativePosition.y;
                plint gZ = iZ + relativePosition.z;

                if ((gX==0 ) | (gY==0 ) | (gZ==0  && Z1>=0) |
                    (gX==X1) | (gY==Y1) | (gZ==Z1 && Z1>=0) ) {
                    topology.get(iX,iY,iZ) = flagConstraintSolid;
                    continue;
                }

                if ((gX==1 ) | (gY==1 ) | (gZ==1  && Z1>=0) |
                    (gX==X1-1) | (gY==Y1-1) | (gZ==Z1-1 && Z1>=0) ) {
                    topology.get(iX,iY,iZ) = flagSolid;
                    continue;
                }

                if (openFlow){
                    topology.get(iX,iY,iZ) = flagBulk;
                } else {
                    topology.get(iX,iY,iZ) = flagSolid;
                }
            }
        }
    }

    if (shape == "d") {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                    // Convert local coordinates to global ones.
                    plint gX = iX + relativePosition.x;
                    plint gY = iY + relativePosition.y;
                    plint gZ = iZ + relativePosition.z;

                    if (gY==Y1 || gY==0) {
                        topology.get(iX,iY,iZ) = flagInlet;
                        continue; }
                    if (gX==X1 || gX==0) {
                        topology.get(iX,iY,iZ) = flagOutlet;
                        continue; }

                    // Convert global coordinates to rotated frame of reference
                    std::vector<float> rotatedCoords = rotate(gX,gY, -1*M_PI/2);
                    plint rX=rotatedCoords[0];
                    plint rY=rotatedCoords[1];

                    // tube connecting in- and outlet
                    if (isInTube(gZ,rX,rY, InRSq,Z1/2,(plint)0, -Y1,Y1*2)) {
                            topology.get(iX,iY,iZ) = flagBulk;
                    }
                }
            }
        }

    }

    if (shape == "Y") {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                    // Convert local coordinates to global ones.
                    plint gX = iX + relativePosition.x;
                    plint gY = iY + relativePosition.y;
                    plint gZ = iZ + relativePosition.z;

                    // circular tube coming from inlet
                    plint Ybifurcation = Y1*2/3;

                    if (heaviside(gY, Ybifurcation)) {
                        if (heaviside(gY, Y1-bufferShiftIn)) {
                            if (isInTube(gZ,gX,gY, InRSq,Z1/2,InC, Ybifurcation,Y1)) {
                                topology.get(iX,iY,iZ) = flagInlet;
                                continue;
                            }
                            if (gX == 0 || gX == X1
                             || gY == 0 || gY == Y1
                             || gZ == 0 || gZ == Z1)
                                continue;

                            if (gY == Y1-bufferShiftIn) {
                                topology.get(iX,iY,iZ) = flagConstraintSolid;
                            } else {
                                topology.get(iX,iY,iZ) = flagSolid;
                            continue;
                            }
                        } else {
                            if (!openFlow) {
                                if (isInTube(gZ,gX,gY, InRSq,Z1/2,InC, (plint)0,Y1)) {
                                    topology.get(iX,iY,iZ) = flagBulk;
                                } else {
                                    topology.get(iX,iY,iZ) = flagSolid;
                                }
                            }
                        }
                    }

                    if (heaviside(gY, Ybifurcation, false)) {
                        if (heaviside(gY, bufferShiftOut, false)) {
                            // left outlet
                            // right outlet
                            if (isInTube(gZ,gX,gY, OutR2Sq,OutRSq, Z1/2,X1/2-OutC, (plint)0,Ybifurcation) ||
                                isInTube(gZ,gX,gY, OutR2Sq,OutRSq, Z1/2,X1/2+OutC, (plint)0,Ybifurcation)) {
                                topology.get(iX,iY,iZ) = flagOutlet;
                                continue;
                            }
                            if (gX == 0 || gX == X1
                             || gY == 0 || gY == Y1
                             || gZ == 0 || gZ == Z1)
                                continue;

                            if (gY == bufferShiftOut) {
                                topology.get(iX,iY,iZ) = flagSolid;
                            } else {
                                topology.get(iX,iY,iZ) = flagConstraintSolid;
                            continue;
                            }
                        } else {
                            if (!openFlow) {
                                if (isInTube(gZ,gX,gY, OutR2Sq,OutRSq, Z1/2,X1/2-OutC, (plint)0,Ybifurcation) ||
                                    isInTube(gZ,gX,gY, OutR2Sq,OutRSq, Z1/2,X1/2+OutC, (plint)0,Ybifurcation)) {
                                        topology.get(iX,iY,iZ) = flagBulk;
                                } else {
                                    topology.get(iX,iY,iZ) = flagSolid;
                                }
                            }
                        }
                    }

                    // horizontal connection
                    if (!openFlow) {
                        if (isInTube(gY,gZ,gX, OutRSq,OutR2Sq, Y1*2/3,Z1/2, X1/2-OutC,X1/2+OutC)) {
                            topology.get(iX,iY,iZ) = flagBulk;
                        }
                    }

                }
            }
        }
    }


    if (shape == "T" || shape == "L" || shape == "Ti") {

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                    // Convert local coordinates to global ones.
                    plint gX = iX + relativePosition.x;
                    plint gY = iY + relativePosition.y;
                    plint gZ = iZ + relativePosition.z;

                    // plane where tubes intersect
                    std::vector<double> normedVec = normalizeVector(-1,1,0);
                    double d = dotProduct(normedVec[0],normedVec[1],normedVec[2], double(InC), double(OutC), 0.);

                    // circular tube coming from inlet
                    plint Ystart;
                    if (shape == "Ti") Ystart = 0;
                    //else              Ystart = OutC;
                    else              Ystart = OutC-OutR;
                    if (isInTube(gZ,gX,gY, InRSq,Z1/2,InC, Ystart,Y1)) {
                        if (gY==Y1 || gY==0) {
                            topology.get(iX,iY,iZ) = flagInlet;
                            continue; }
                        // upper inlet buffer
                        if (heaviside(gY, Y1-bufferShiftIn)) {
                            topology.get(iX,iY,iZ) = flagBuffer;
                            continue; }
                        if (shape == "Ti") {
                            // lower inlet buffer
                            if (heaviside(gY, 0+bufferShiftIn, false)) {
                                topology.get(iX,iY,iZ) = flagBuffer;
                                continue; }
                        }
                        if (shape=="T" || isBelowPlane(double(gX),double(gY),double(gZ), normedVec[0],normedVec[1],normedVec[2], d) == false ) {
                            topology.get(iX,iY,iZ) = flagBulk;
                            continue; }
                    } else {
                        if (heaviside(gY, Y1-bufferShiftIn)) {
                            topology.get(iX,iY,iZ) = flagConstraintSolid;
                            continue; }
                        if (shape == "Ti") {
                            // lower inlet buffer
                            if (heaviside(gY, 0+bufferShiftIn, false)) {
                                topology.get(iX,iY,iZ) = flagConstraintSolid;
                                continue; }
                        }
                    }

                    // circular tube coming from outlet
                    plint Xstart;
                    if (shape == "T") Xstart = 0;
                    //else              Xstart = InC;
                    else              Xstart = InC-InR;
                    if (isInTube(gY,gZ,gX, OutRSq,OutR2Sq, OutC,Z1/2, Xstart,X1)) {
                        if (gX==X1 || gX==0) {
                            topology.get(iX,iY,iZ) = flagOutlet;
                            continue; }
                        // right outlet buffer
                        if (heaviside(gX, X1-bufferShiftOut)) {
                            topology.get(iX,iY,iZ) = flagBuffer;
                            continue; }
                        if (shape == "T") {
                            // left outlet buffer
                            if (heaviside(gX, 0+bufferShiftOut, false)) {
                                topology.get(iX,iY,iZ) = flagBuffer;
                                continue; }
                        }
                        if (shape=="T" || isBelowPlane(double(gX),double(gY),double(gZ), normedVec[0],normedVec[1],normedVec[2], d) == true ) {
                            topology.get(iX,iY,iZ) = flagBulk;
                            continue; }
                    } else {
                        if (heaviside(gX, X1-bufferShiftOut)) {
                            topology.get(iX,iY,iZ) = flagConstraintSolid;
                            continue; }
                        if (shape == "T") {
                            if (heaviside(gX, 0+bufferShiftOut, false)) {
                                topology.get(iX,iY,iZ) = flagConstraintSolid;
                                continue; }
                        }
                    }
                }
            }
        }
    }

    if (shape == "-") {

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                    // Convert local coordinates to global ones.
                    plint gX = iX + relativePosition.x;
                    plint gY = iY + relativePosition.y;
                    plint gZ = iZ + relativePosition.z;

                    // circular tube coming from inlet
                    if (isInTube(gY,gZ,gX, InRSq,InC,Z1/2, (plint)(0),X1/8)) {
                        if (gX==0) {
                            topology.get(iX,iY,iZ) = flagInlet;
                            continue; }
                        // inlet buffer
                        if (heaviside(gX, 0+bufferShiftIn, false)) {
                            topology.get(iX,iY,iZ) = flagBuffer;
                            continue; }
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    } else {
                        if (heaviside(gX, 0+bufferShiftIn, false)) {
                            topology.get(iX,iY,iZ) = flagConstraintSolid;
                            continue; }
                    }
                    // circular tube coming from outlet
                    if (isInTube(gY,gZ,gX, InRSq,InC,Z1/2, X1/4,X1)) {
                        if (gX==X1) {
                            topology.get(iX,iY,iZ) = flagOutlet;
                            continue; }
                        // right outlet buffer
                        if (heaviside(gX, X1-bufferShiftOut)) {
                            topology.get(iX,iY,iZ) = flagBuffer;
                            continue; }
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    } else {
                        if (heaviside(gX, X1-bufferShiftOut)) {
                            topology.get(iX,iY,iZ) = flagConstraintSolid;
                            continue; }
                    }
                    // circular tube shifted || (shiftedTube)
                    if (isInTube(gY,gZ,gX, OutRSq,OutC,Z1/2, X1/8,X1/4)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // circular tube going up near inlet (diversionU)
                    if (isInTube(gZ,gX,gY, ((InR>OutR)? InRSq : OutRSq),Z1/2,X1/8, InC,OutC)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // circular tube going down near outlet (diversionD)
                    if (isInTube(gZ,gX,gY, ((InR>OutR)? InRSq : OutRSq),Z1/2,X1/4, InC,OutC)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // spherical part where cylinders meet (inlet/diversionU)
                    if (isInSphere(gX,gY,gZ, ((InR>OutR)? InRSq : OutRSq), X1/8,InC,Z1/2)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // spherical part where cylinders meet (diversionU/shiftedTube)
                    if (isInSphere(gX,gY,gZ, ((InR>OutR)? InRSq : OutRSq), X1/8,OutC,Z1/2)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // spherical part where cylinders meet (shiftedTube/diversionD)
                    if (isInSphere(gX,gY,gZ, ((InR>OutR)? InRSq : OutRSq), X1/4,OutC,Z1/2)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // spherical part where cylinders meet (diversionU/outlet)
                    if (isInSphere(gX,gY,gZ, ((InR>OutR)? InRSq : OutRSq), X1/4,InC,Z1/2)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }

                }
            }
        }
    }

    if (shape == "U") {

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                    // Convert local coordinates to global ones.
                    plint gX = iX + relativePosition.x;
                    plint gY = iY + relativePosition.y;
                    plint gZ = iZ + relativePosition.z;

                    // circular tube coming from inlet
                    if (isInTube(gZ,gX,gY, InR2Sq,InRSq, Z1/2,InC, Y1/4-cross,Y1)) {
                        if (gY==Y1) {
                            topology.get(iX,iY,iZ) = flagInlet;
                            continue; }
                        // inlet buffer
                        if (heaviside(gY, Y1-bufferShiftIn)) {
                            topology.get(iX,iY,iZ) = flagBuffer;
                            continue; }
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // circular tube coming from outlet
                    if (isInTube(gZ,gX,gY, OutR2Sq,OutRSq, Z1/2,OutC, Y1/4-cross,Y1)) {
                        if (gY==Y1) {
                            topology.get(iX,iY,iZ) = flagOutlet;
                            continue; }
                        // outlet buffer
                        if (heaviside(gY, Y1-bufferShiftOut)) {
                            topology.get(iX,iY,iZ) = flagBuffer;
                            continue; }
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    if (!isInTube(gZ,gX,gY, OutR2Sq,OutRSq, Z1/2,OutC, Y1/4,Y1) &&
                        !isInTube(gZ,gX,gY, InR2Sq,InRSq, Z1/2,InC, Y1/4,Y1)) {
                        if (heaviside(gY, Y1-bufferShiftIn)) {
                            topology.get(iX,iY,iZ) = flagConstraintSolid;
                            continue; }
                    }
                    // circular tube X-parallel
                    if (isInTube(gY,gZ,gX, InRSq,InR2Sq, Y1/4,Z1/2, InC-cross,OutC+cross)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    /*
                    // spherical part where inlet meets X-parallel tube
                    if (isInSphere(gX,gY,gZ, ((InR>OutR)? InRSq : OutRSq), InC,Y1/4,Z1/2)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    // spherical part where X-parallel tube meets outlet
                    if (isInSphere(gX,gY,gZ, ((InR>OutR)? InRSq : OutRSq), OutC,Y1/4,Z1/2)) {
                        topology.get(iX,iY,iZ) = flagBulk;
                        continue;
                    }
                    */
                }
            }
        }

    }

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                int& cell = topology.get(iX,iY,iZ);
                if (cell == flagConstraintSolid) {
                    // find 'constraintBB' nodes
                    if (isAnyDirectNeighbourValue3D(topology, iX,iY,iZ, flagBuffer)) {
                        cell = flagConstraintBB;
                        continue;
                    }
                }
            }
        }
    }
}

/* ******** PPM2TopologyFunctional3D ************************************** */
template<typename T>
PPM2TopologyFunctional3D<T>::PPM2TopologyFunctional3D(std::string PPMname_)
    : PPMname(PPMname_)
{
    iPPM = easyppm_create(1, 1, IMAGETYPE_PPM);
    easyppm_read(&iPPM, PPMname.c_str());
    easyppm_invert_y(&iPPM);
}

template<typename T>
void PPM2TopologyFunctional3D<T>::process (
        Box3D domain,
        ScalarField3D<T>& topology)
{
    Dot3D relativePosition = topology.getLocation();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // Convert local coordinates to global ones.
            plint i = iX + relativePosition.x;
            plint j = iY + relativePosition.y;
//            pcout << "topology(" << iX << "," << iY <<") = " << topology.get(iX,iY);
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
                topology.get(iX,iY,iZ) = (int)(easyppm_get(&iPPM, i,j).r);
//            pcout << " -> " << topology.get(iX,iY)
//                  << " (" << (int)(easyppm_get(&iPPM, i,j).r) << ")" << endl;
        }
    }
}


/* ******** UpdateTopologyFunctional3D ************************************** */
template<typename T>
void UpdateTopologyFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& topology )
{

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                int& cell = topology.get(iX,iY,iZ);

                if (cell == flagSolid) {
                    // find 'BB' nodes
                    // also look for diagonally neighbouring fluid nodes!
                    if (   (isAnyDirectNeighbourValue3D(topology, iX,iY,iZ, flagWet)
                         || isAnyDirectNeighbourValue3D(topology, iX,iY,iZ, flagBulk))
                        && (isAnyDirectNeighbourValue3D(topology, iX,iY,iZ, flagSolid)
                         || isAnyDirectNeighbourValue3D(topology, iX,iY,iZ, flagConstraintSolid))) {
                        cell = flagSedimenting;
                        continue;
                    }
                }

                if (cell == flagBB) {
                    // in case a bbNode ends up surrounded by solids
                    if (!isAnyDirectNeighbourValue3D(topology, iX,iY,iZ, flagWet)
                        && !isAnyDirectNeighbourValue3D(topology, iX,iY,iZ, flagBulk) ) {
                        cell = flagSolid;
                        continue;
                    }
                }

                if (cell == flagBulk) {
                    // find 'wet' nodes
                    if (isAnyNearestNeighbourValue3D(topology, iX,iY,iZ, flagBB)
                        || isAnyNearestNeighbourValue3D(topology, iX,iY,iZ, flagSedimenting)
                        || isAnyNearestNeighbourValue3D(topology, iX,iY,iZ, flagSolid)) {
                        cell = flagWet;
                        continue;
                    }
                }

                if (cell == flagWet) {
                // in case a wetNode loses its neighbouring BBs
                    if (!isAnyNearestNeighbourValue3D(topology, iX,iY,iZ, flagBB)
                        && !isAnyNearestNeighbourValue3D(topology, iX,iY,iZ, flagSedimenting)
                        && !isAnyNearestNeighbourValue3D(topology, iX,iY,iZ, flagSolid)) {
                        cell = flagBulk;
                        continue;
                    }
                }

                /*
                if (cell == flagBulk) {
                    // in case a bulkNode ends up surrounded by BBs
                    if (isAllNearestNeighbourValue3D(topology, iX,iY,iZ, flagBB)
                        ||isAllNearestNeighbourValue3D(topology, iX,iY,iZ, flagSedimenting)
                        ||isAllNearestNeighbourValue3D(topology, iX,iY,iZ, flagSolid)) {
                        cell = flagWet;
                        continue;
                    }
                }
                */

            }
        }
    }

}

/* ******** FreezeTopologyFunctional3D ************************************** */
template<typename T>
void FreezeTopologyFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& topology )
{

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                int& cell = topology.get(iX,iY,iZ);

                if (cell == flagBulk) {
                    cell = flagBuffer;
                }
                if (cell == flagBB) {
                    cell = flagConstraintBB;
                }
            }
        }
    }

}

/* ******** SedimentFunctional3D ********************************************** */
template<typename T1, typename T2>
SedimentFunctional3D<T1,T2>::SedimentFunctional3D(T2 sediThreshold_)
    : sediThreshold(sediThreshold_)
{ }

template<typename T1, typename T2>
void SedimentFunctional3D<T1,T2>::process (
        Box3D domain, ScalarField3D<T1>& topology, ScalarField3D<T2>& shearStressField )
{
    Dot3D offset = computeRelativeDisplacement(topology, shearStressField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                // lessEqual sediThreshold && flagWet?
                // -> let them sediment
                //      mark them as flagSedimenting
                //      update dynamics to BB
                //      and then mark them as flagBB

                if (topology.get(iX,iY,iZ)==flagWet) {
                    if (shearStressField.get(iX+offset.x,iY+offset.y,iZ+offset.z)<=sediThreshold) {
                        topology.get(iX,iY,iZ) = flagSedimenting;
                    }
                }
            }
        }
    }

}

/* ******** ErodeFunctional3D ********************************************** */
template<typename T1, typename T2>
void ErodeFunctional3D<T1,T2>::process (
        Box3D domain, ScalarField3D<T1>& topology, ScalarField3D<T2>& shearStressField )
{
    Box3D extendedDomain(domain.enlarge(1));

    Dot3D offset = computeRelativeDisplacement(topology, shearStressField);
    for (plint iX=extendedDomain.x0; iX<=extendedDomain.x1; ++iX) {
        for (plint iY=extendedDomain.y0; iY<=extendedDomain.y1; ++iY) {
            for (plint iZ=extendedDomain.z0; iZ<=extendedDomain.z1; ++iZ) {

                // greaterEqual erosThreshold && flagWet?
                // -> let them erode
                //      check neighbours if they are flagBB
                //      -> mark them as flagEroding (and make them fluid nodes after functional finished)

                if (topology.get(iX,iY,iZ)==flagWet) {
                    if (shearStressField.get(iX+offset.x,iY+offset.y,iZ+offset.z)>=erosThreshold) {

                        for(int x=-1;x<=1;x++)
                            for(int y=-1;y<=1;y++)
                                for(int z=-1;z<=1;z++) {
                                    if (x==0 && y==0 && z==0) continue;
                                    if (sign(x)+sign(y)+sign(z) > 1) continue;
                                    if (topology.get(iX+x,iY+y,iZ+z) == flagBB ||
                                        topology.get(iX+x,iY+y,iZ+z) == flagWet)
                                        topology.get(iX+x,iY+y,iZ+z) = flagEroding;
                                }
                        topology.get(iX,iY,iZ) = flagBulk;

                    }
                }

            }
        }
    }

}

namespace plb {

template<typename T>
TopoWriter<T>::TopoWriter(plint valueRange_)
    : valueRange(valueRange_)
{ }

template<typename T>
void TopoWriter<T>::writePpm (
        std::string const& fName, bool colored,
        MultiScalarField2D<T>& field) const
{
    global::profiler().start("io");
    ScalarField2D<T> localField(field.getNx(), field.getNy());
    copySerializedBlock(field, localField);
    writePpmImplementation(fName, colored, localField);
    global::profiler().stop("io");
}

template<typename T>
void TopoWriter<T>::writePpm (
        std::string const& fName, bool colored,
        MultiScalarField3D<T>& field, int dir) const
{
    global::profiler().start("io");
    ScalarField2D<T> localField = extractMiddleLayer(field, dir);
    writePpmImplementation(fName, colored, localField);
    global::profiler().stop("io");
}

template<typename T>
void TopoWriter<T>::writePpmImplementation (
        std::string const& fName, bool colored,
        ScalarField2D<T>& localField) const
{
    if (global::mpi().isMainProcessor()) {
        T maxVal = computeMax(localField);
        T minVal = computeMin(localField);
        if (maxVal > valueRange || minVal < 0.) {
            pcout << "ScalarField has values that exceed valueRange!" << endl
                  << "maxVal = " << maxVal << ", minVal = " << minVal << endl
                  << "No ppm written." << endl;
            return;
        }
        std::string fullName;
        if (colored)
            fullName = global::directories().getImageOutDir() + fName+".ppm";
        else
            fullName = global::directories().getImageOutDir() + fName+".pgm";
        std::ofstream fout(fullName.c_str());
        if (colored) fout << "P3\n";
        else         fout << "P2\n";
        fout << localField.getNx() << " " << localField.getNy() << "\n";
        fout << (valueRange-1) << "\n";

        for (plint iY=localField.getNy()-1; iY>=0; --iY) {
            for (plint iX=0; iX<localField.getNx(); ++iX) {
                // create variables to hold rgb value for PPM export:
                double red, green, blue;
                // red always represents a flag value of topology
                red = (double) (localField.get(iX,iY));

                switch( (int)(red) ) {
                    // for these flags, colors are nice to spot them easily
                    case flagSedimenting:     green= 255.; blue=red= 50.; break;
                    case flagEroding:         green= 255.; blue=red= 0.;  break;
                    case flagBuffer:          green= 0.;   blue=   0.;    break;
                    case flagWet:             green= 200.; blue=   0.;    break;
                    case flagConstraintSolid: green= 0.;   blue= 255.;    break;
                    case flagConstraintBB:    green= 50.;  blue= 255.;    break;
                    // BB and Bulk are not changed, rgb is all set to same value
                    default: green=blue=red;
                }

                if (colored) {
                    fout << (int) (red)   << " "
                         << (int) (green) << " "
                         << (int) (blue)  << "\n";
                } else {
                    fout << (int) (red)   << "\n"; // greyValue
                }
            }
        }
        fout.close();
    }
}

template<typename T>
void TopoWriter<T>::writeBinary (
        std::string const& fName,
        MultiScalarField3D<T>& field) const
{
    global::profiler().start("io");
    plint nx=field.getNx();
    plint ny=field.getNy();
    plint nz=field.getNz();
    ScalarField3D<T> localField(nx,ny,nz);
    copySerializedBlock(field, localField);
    writeBinaryImplementation(fName, localField);
    global::profiler().stop("io");
}

template<typename T>
void TopoWriter<T>::writeBinaryImplementation (
        std::string const& fName,
        ScalarField3D<T>& localField) const
{
    if (global::mpi().isMainProcessor()) {
        T maxVal = computeMax(localField);
        T minVal = computeMin(localField);
        if (maxVal > valueRange || minVal < 0.) {
            pcout << "ScalarField has values that exceed valueRange!" << endl
                  << "maxVal = " << maxVal << ", minVal = " << minVal << endl
                  << "No file written." << endl;
            return;
        }
        std::string fullName;
        fullName = global::directories().getImageOutDir() + fName+".binary";

        std::ofstream fout(fullName.c_str(), std::ios::out | std::ios::binary);

        for (plint iZ=0; iZ<localField.getNz(); ++iZ) {
            //for (plint iY=localField.getNy()-1; iY>=0; --iY) {
            for (plint iY=0; iY<localField.getNy(); ++iY) {
                for (plint iX=0; iX<localField.getNx(); ++iX) {
                    char value = localField.get(iX,iY,iZ);
                    fout.write( (char*) &value, sizeof(value) );

                    //T value = localField.get(iX,iY,iZ);
                    //fout.write( (char*) &value, sizeof(T) );
                }
            }
        }
        fout.close();
    }
}
}

#endif
