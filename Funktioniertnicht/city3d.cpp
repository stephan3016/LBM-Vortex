/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#include <olb3D.h>
#include <olb3D.hh>

using namespace olb;

using T = float;
using D = double;
using DESCRIPTOR = descriptors::D3Q19<descriptors::POROSITY>;
using BulkDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>
>;

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& sGeometry,
                     IndicatorF3D<T>& volume)
{
  sGeometry.rename(0,2);
  sGeometry.rename(2,1,{1,1,1});

    //Floor: mat_num 3
  {
    auto origin = volume.getMin() - 1*converter.getPhysDeltaX();
    auto extent = volume.getMax() + 1*converter.getPhysDeltaX();
    extent[2] = 0.2*volume.getMax()[2];
    IndicatorCuboid3D<T> floor(extent, origin);
    sGeometry.rename(2,3,floor);
  }

    // Inlet: mat_num 4
  {
      auto boundaryVector = volume.getMax()+1*converter.getPhysDeltaX();
      boundaryVector[0] = 1*converter.getPhysDeltaX();

      auto boundaryOrigin = volume.getMin() - 1*converter.getPhysDeltaX();
      IndicatorCuboid3D<T> VortexBoundary(boundaryVector, boundaryOrigin);
      sGeometry.rename(2,4,1, VortexBoundary);
  }

    // outlet: mat_num 5
    {
        auto TopLayerVector = volume.getMax() + 1*converter.getPhysDeltaX();

        auto TopLayerOrigin = volume.getMin() - 1*converter.getPhysDeltaX();
        IndicatorCuboid3D<T> TopLayerBoundary(TopLayerVector, TopLayerOrigin);
        sGeometry.rename(2,5,1, TopLayerBoundary);
    }

    // Hence, the top+side layer is mat_num 2
  sGeometry.clean();
  sGeometry.innerClean();
  sGeometry.checkForErrors();
  sGeometry.print();
}

void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,3>& sGeometry,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    IndicatorF3D<T>& volume)
{
  const T omega = converter.getLatticeRelaxationFrequency();

  sLattice.defineDynamics<NoDynamics>(sGeometry, 0);
  sLattice.defineDynamics<BulkDynamics>(sGeometry, 1);
  sLattice.defineDynamics<BounceBack>(sGeometry, 2);
  sLattice.defineDynamics<BounceBack>(sGeometry, 3);
  sLattice.defineDynamics<BulkDynamics>(sGeometry, 4);
  sLattice.defineDynamics<BulkDynamics>(sGeometry, 5);
  sLattice.defineDynamics<BulkDynamics>(sGeometry, 6);



  {
    setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, sGeometry, 4);
    setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, sGeometry, 5);
    setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, sGeometry, 6);
    AnalyticalConst3D<T,T> rhoF(1);
    AnalyticalConst3D<T,T> uF(0, 0, 0);
    sLattice.defineRhoU(sGeometry.getMaterialIndicator(2), rhoF, uF);
  }

  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.setParameter<collision::LES::Smagorinsky>(0.1);

  {
    auto bulkIndicator = sGeometry.getMaterialIndicator({1,3});
    AnalyticalConst3D<T,T> rhoF( T( 1 ) );
    AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) );
    sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);
    sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  }

  AnalyticalConst3D<T,T> initialPorosityF(1);
  sLattice.defineField<descriptors::POROSITY>(sGeometry.getMaterialIndicator({0,1,2,3,4,5,6}), initialPorosityF);

  AnalyticalConst3D<T,T> solidPorosityF(0);
  SuperIndicatorFfromIndicatorF3D<T> discreteVolume(
    volume, sGeometry);
  sLattice.defineField<descriptors::POROSITY>(discreteVolume, solidPorosityF);

  sLattice.initialize();
}

void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
                       SuperGeometry<T,3>& sGeometry,
                       SuperLattice<T,DESCRIPTOR>& sLattice,
                       std::size_t iT,  IndicatorF3D<T>& volume)
{
  OstreamManager clout(std::cout, "boundary");

  const auto maxStartT =  converter.getLatticeTime(60);
  const auto startIterT = converter.getLatticeTime(0.1);

  if (iT < maxStartT && iT % startIterT == 0) {

    PolynomialStartScale<T,std::size_t> scale(maxStartT, 1);
    T frac{};
    scale(&frac, &iT);

    // paramters for exponential function:
    std::vector<T> alpha(3,0);
    std::vector<T> beta(3,0);
    std::vector<T> betaVel(3,0);
    betaVel[0] = 0.004*frac;
    beta[0] = 100;
    alpha[0] = 4;



    //Implementation of the vortex method:

// -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*
     //1: Arguments for the vortex method
if(0){//VortexMethode ein/aus
    // a) Construct an IndicatorF3D in the form of an IndicatorCuboid3D for the inlet:
      auto boundaryVector = volume.getMax() -volume.getMin()+ 1*converter.getPhysDeltaX();
      boundaryVector[0] = 1*converter.getPhysDeltaX();
      auto boundaryOrigin = volume.getMin() - 1*converter.getPhysDeltaX();
      //std::cout<<"XXXXXXXXXXXXXXXXXXXXorigin"<<boundaryOrigin<<"  Vektor "<<boundaryVector<<std::endl;
      boundaryOrigin[2] = 0.2*volume.getMax()[2];
      IndicatorCuboid3D<T> VortexBoundaryIndicator(boundaryVector, boundaryOrigin);

      // b) declare arguments of the type FunctorPtr<...>
      auto VortexBoundary = FunctorPtr<IndicatorF3D<T>>(new IndicatorCuboid3D<T>(boundaryVector, boundaryOrigin));
      auto InletLatticeI = FunctorPtr<SuperIndicatorF3D<T>>(new SuperIndicatorFfromIndicatorF3D<T>(VortexBoundaryIndicator, sGeometry));

      // c) declare the direction vector of type T (=float)
      Vector<T,3> axisDirection(T(1),T(0),T(0));



    //2: Declare the Vortex method as a class by calling its constructor
    // In this case, the class will be called "VortexMethodTurbulentVelocityBoundary"
    VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR> VortexMethodTurbulentVelocityBoundary(std::move(InletLatticeI), std::move(VortexBoundary), converter, sLattice,42,100, (D) 0.16, (D) 0.9, axisDirection);



    //3. Implementing the desired velocity profile for the vortex boundary
    // a) The vortex method requires the profile to be given not as AnalyticalF3D but as a shared_ptr pointing to the actual AnalyticalF3D that indicates the profile
    auto uF = std::shared_ptr<AnalyticalF3D<T,T>>(new ExponentialProfile3D<T>(sGeometry,4, alpha, beta, betaVel));

    // b) Now we call the functions inside the class we created above
    VortexMethodTurbulentVelocityBoundary.setProfile(uF);       // to set the desired v-profile



    //4. We apply all changes to the boundary
    VortexMethodTurbulentVelocityBoundary.apply(iT);            // to apply our changes



    // 5. VortexBoundary condition was applied :)

// -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*
}
else{
    ExponentialProfile3D<T> ExpProfile4(sGeometry,4, alpha, beta, betaVel);
    clout << iT << " " << frac << std::endl;
    sLattice.defineU(sGeometry, 4, ExpProfile4);
}

    // Define the exponential profile in z-direction. The resulting object is of type AnalyticalF3D
    // This here is only for the top boundary layer, where no Vortex method is implemented so far --> is it necessary to implemented there?
    ExponentialProfile3D<T> ExpProfile2(sGeometry,5, alpha, beta, betaVel);
    sLattice.defineU(sGeometry, 5, ExpProfile2);
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(std::size_t iT,
                SuperLattice<T,DESCRIPTOR>& sLattice,
                SuperGeometry<T,3>& sGeometry,
                const UnitConverter<T,DESCRIPTOR>& converter,
                util::Timer<T>& timer)
{
  SuperVTMwriter3D<T> vtmWriter("city3d");

  if (iT == 0) {
    std::remove("output.dat");
    SuperLatticeGeometry3D<T,DESCRIPTOR> geometryF(sLattice, sGeometry);
    SuperLatticeCuboid3D<T,DESCRIPTOR> cuboidF(sLattice);
    SuperLatticeRank3D<T,DESCRIPTOR> rankF(sLattice);
    SuperLatticeField3D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(sLattice);

    vtmWriter.write(geometryF);
    vtmWriter.write(cuboidF);
    vtmWriter.write(rankF);
    vtmWriter.write(porosityF);
    vtmWriter.createMasterFile();
  }
  std::ofstream fout("output.dat", std::ios::app);					//Ausgabe Stats in datei
 fout<<"Average Rho+"<<sLattice.getStatistics().getAverageRho()<<"+Average Energy+"<<sLattice.getStatistics().getAverageEnergy()<<std::endl;

  if (iT % converter.getLatticeTime(1) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  if (iT % converter.getLatticeTime(1) == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("city3d");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      //SuperLatticeDensity3D density(sLattice, converter);
      vtkWriter.addFunctor(velocity);
      //vtkWriter.addFunctor(density);
      task(vtkWriter, iT);
    });
  }
}

int main(int argc, char* argv[]) {
  olbInit(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  CLIreader args(argc, argv);

  const UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
    int {1},     // resolution: number of voxels per charPhysL
    (T)   0.1,   // lattice velocity
    (T)   1.0,   // charPhysLength: reference length of simulation geometry
    (T)  50.0,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   10e-5, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  // Available at https://bwsyncandshare.kit.edu/s/zHRZSHnRCwAw5QA
  STLreader<T> volume2("Mensa_klein.stl", converter.getCharPhysLength(), 1, 1);
  int einlauflaenge=100;
  auto extent = volume2.getMax() - volume2.getMin();
  extent[0]+=einlauflaenge*converter.getPhysDeltaX();
  extent[2] *= 2;

  
  
  //extent[0] += 100*converter.getPhysDeltaX();//laengere Box fuer Einlauf
  auto Position=volume2.getMin();
  Position[0]-=einlauflaenge*converter.getPhysDeltaX();
  //IndicatorCuboid3D<T> bounding(extent, Position);
  
  IndicatorCuboid3D<T> volume(extent, Position);
  CuboidGeometry3D<T> cGeometry(volume, converter.getPhysDeltaX(), 1*singleton::mpi().getSize());
  //cGeometry.setPeriodicity(true, true, false);

  BlockLoadBalancer<T> loadBalancer(cGeometry);
  SuperGeometry<T,3> sGeometry(cGeometry, loadBalancer);

  prepareGeometry(converter, sGeometry, volume);

  SuperLattice<T,DESCRIPTOR> sLattice(sGeometry);

  prepareLattice(sLattice, sGeometry, converter, volume);

  util::Timer<T> timer(converter.getLatticeTime(1e6), sGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT <= converter.getLatticeTime(1e6); ++iT) {
    setBoundaryValues(converter, sGeometry, sLattice, iT, volume);
    sLattice.collideAndStream();
    getResults(iT, sLattice, sGeometry, converter, timer);
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();

  return 0;
}
