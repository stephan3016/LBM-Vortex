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
using DESCRIPTOR = descriptors::D3Q19<descriptors::POROSITY>;
using BulkDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>
>;


void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& sGeometry,
                     IndicatorF3D<T>& cavity)
{
    sGeometry.rename(0,2);
    sGeometry.rename(2,1,{1,1,1});


    auto origin = cavity.getMin()- 1*converter.getPhysDeltaX();
    auto extent = cavity.getMax()+ 1*converter.getPhysDeltaX();
    extent[0] = 1*converter.getPhysDeltaX();
    IndicatorCuboid3D<T> inflow(extent,origin);
    sGeometry.rename(2,3,1, inflow);

    origin[0] = sGeometry.getStatistics().getMaxPhysR( 2 )[0]-converter.getConversionFactorLength();
    IndicatorCuboid3D<T> outflow(extent,origin);
    sGeometry.rename(2,4,1, outflow);


    sGeometry.clean();
    sGeometry.innerClean();
    sGeometry.checkForErrors();
    sGeometry.print();


}


void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,3>& sGeometry,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    IndicatorF3D<T>& cavity)
{
    const T omega = converter.getLatticeRelaxationFrequency();

    sLattice.defineDynamics<NoDynamics>(sGeometry, 0);
    sLattice.defineDynamics<BulkDynamics>(sGeometry, 1);
    sLattice.defineDynamics<BulkDynamics>(sGeometry, 3);
    sLattice.defineDynamics<BulkDynamics>(sGeometry, 4);
    sLattice.defineDynamics<BulkDynamics>(sGeometry, 2);

    setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, sGeometry, 2);       // Haftbedingung
    AnalyticalConst3D<T,T> rhoF(1);
    AnalyticalConst3D<T,T> uF(0, 0, 0);
    sLattice.defineRhoU(sGeometry.getMaterialIndicator(2), rhoF, uF);

    sLattice.setParameter<descriptors::OMEGA>(omega);
    sLattice.setParameter<collision::LES::Smagorinsky>(0.1);


    auto bulkIndicator = sGeometry.getMaterialIndicator({1,3,4});
    AnalyticalConst3D<T,T> rhoF2( T( 1 ) );
    AnalyticalConst3D<T,T> uF2( T( 0 ), T( 0 ), T( 0 ) );
    sLattice.iniEquilibrium(bulkIndicator, rhoF2, uF2);
    sLattice.defineRhoU(bulkIndicator, rhoF2, uF2);

    AnalyticalConst3D<T,T> solidPorosityF(0);
  SuperIndicatorFfromIndicatorF3D<T> discreteVolume(cavity, sGeometry);
  sLattice.defineField<descriptors::POROSITY>(discreteVolume, solidPorosityF);

  sLattice.initialize();
}


void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
                       SuperGeometry<T,3>& sGeometry,
                       SuperLattice<T,DESCRIPTOR>& sLattice,
                       std::size_t iT)
{
  OstreamManager clout(std::cout, "boundary");

  const auto maxStartT =  converter.getLatticeTime(60);
  const auto startIterT = converter.getLatticeTime(0.1);

  if (iT < maxStartT && iT % startIterT == 0) {
    auto uF = std::shared_ptr<AnalyticalF3D<T,T>>(new AnalyticalConst3D<T,T>(converter.getCharLatticeVelocity(), 0, 0));

    PolynomialStartScale<T,std::size_t> scale(maxStartT, 1);
    T frac{};
    scale(&frac, &iT);

    clout << iT << " " << frac << std::endl;

    sLattice.defineU(sGeometry, 3, *(frac * uF));

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
  SuperLatticeGeometry3D<T, DESCRIPTOR> materials( sLattice, sGeometry );
  vtmWriter.addFunctor( materials );

  if (iT == 0) {
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

  if (iT % converter.getLatticeTime(1) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  if (iT % converter.getLatticeTime(60) == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("city3d");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      vtkWriter.addFunctor(velocity);
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
    (T)  10.0,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   10e-5, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  converter.print();



  //Creating the simulation box:
    Vector<T,3> extent(5,5,5);
    Vector<T,3> origin(0,0,0);
    IndicatorCuboid3D<T> cavity(extent,origin);
    CuboidGeometry3D<T> cGeometry(cavity, converter.getPhysDeltaX(), singleton::mpi().getSize());

    BlockLoadBalancer<T> loadBalancer(cGeometry);
    SuperGeometry<T,3> sGeometry(cGeometry, loadBalancer);

    prepareGeometry(converter, sGeometry, cavity);

    SuperLattice<T,DESCRIPTOR> sLattice(sGeometry);

    prepareLattice(sLattice, sGeometry, converter, cavity);

    util::Timer<T> timer(converter.getLatticeTime(1e6), sGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT <= converter.getLatticeTime(1e6); ++iT) {
    setBoundaryValues(converter, sGeometry, sLattice, iT);
    sLattice.collideAndStream();
    getResults(iT, sLattice, sGeometry, converter, timer);
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();

  return 0;
}


