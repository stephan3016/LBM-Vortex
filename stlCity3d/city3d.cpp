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
                     IndicatorF3D<T>& volume)
{
  sGeometry.rename(0,2);
  sGeometry.rename(2,1,{1,1,1});

  {
    auto origin = volume.getMin() - 1*converter.getPhysDeltaX();
    auto extent = volume.getMax() + 1*converter.getPhysDeltaX();
    extent[2] = 0.2*volume.getMax()[2];
    IndicatorCuboid3D<T> floor(extent, origin);
    sGeometry.rename(2,3,floor);
  }

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
  sLattice.defineDynamics<BounceBack>(sGeometry, 3);

  {
    setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, sGeometry, 2);
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
  sLattice.defineField<descriptors::POROSITY>(sGeometry.getMaterialIndicator({0,1,2,3}), initialPorosityF);

  AnalyticalConst3D<T,T> solidPorosityF(0);
  SuperIndicatorFfromIndicatorF3D<T> discreteVolume(
    volume, sGeometry);
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

    sLattice.defineU(sGeometry, 2, *(frac * uF));

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

  // Available at https://bwsyncandshare.kit.edu/s/zHRZSHnRCwAw5QA
  STLreader<T> volume("Mensa_klein.stl", converter.getCharPhysLength(), 1, 1);
  auto extent = volume.getMax() - volume.getMin();
  extent[2] *= 2;
  IndicatorCuboid3D<T> bounding(extent, volume.getMin());
  CuboidGeometry3D<T> cGeometry(bounding, converter.getPhysDeltaX(), 2*singleton::mpi().getSize());
  //cGeometry.setPeriodicity(true, true, false);

  BlockLoadBalancer<T> loadBalancer(cGeometry);
  SuperGeometry<T,3> sGeometry(cGeometry, loadBalancer);

  prepareGeometry(converter, sGeometry, bounding);

  SuperLattice<T,DESCRIPTOR> sLattice(sGeometry);

  prepareLattice(sLattice, sGeometry, converter, volume);

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
