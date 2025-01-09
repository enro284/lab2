#include <array>
#include <iostream>

#include "particle.hpp"

int Particle::numParticleType_{0};
const ParticleType *Particle::particleTypeArray_[maxNumParticleType_];

int main() {
  ParticleType particleType{"particella", 10., 5};
  particleType.Print();
  std::cout << particleType.GetName() << ' ' << particleType.GetMass() << ' '
            << particleType.GetCharge() << '\n';

  ResonanceType resonanceType{"particella con risonanza", 13., 4, 100.};
  resonanceType.Print();
  std::cout << resonanceType.GetName() << ' ' << resonanceType.GetMass() << ' '
            << resonanceType.GetCharge() << ' ' << resonanceType.GetWidth()
            << '\n';

  ParticleType *test[2];
  test[0] = &particleType;
  test[1] = &resonanceType;

  for (int i = 0; i < 2; i++) {
    test[i]->Print();
  }

  Particle::AddParticleType("pi", 12., 3.);
  Particle::AddParticleType("e-", 7., 14.);
  Particle::AddParticleType("mu", 90., 83.);
  Particle particle1{"pi", 15., 18., 32.};
  Particle particle2{"pi", 33., 27., 48};
  Particle::PrintParticleTypeArray();

  std::cout << "energia totale = " << particle1.EnergyTot() << '\n';
  std::cout << "massa invariante = " << particle1.InvMass(particle2) << '\n';
  particle1.SetP(3., 21., 42);
  std::cout << "energia totale = " << particle1.EnergyTot() << '\n';
  std::cout << "massa invariante = " << particle1.InvMass(particle2) << '\n';
}