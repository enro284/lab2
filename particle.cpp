#include "particle.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

const char *ParticleType::GetName() const { return name_; }

const double ParticleType::GetMass() const { return mass_; }

const int ParticleType::GetCharge() const { return charge_; }

void ParticleType::Print() const {
  std::cout << "Name: " << name_ << ", Mass: " << mass_
            << ", Charge: " << charge_ << '\n';
}
void ResonanceType::Print() const {
  ParticleType::Print();
  std::cout << ", Width: " << width_ << '\n';
}

double ResonanceType::GetWidth() const { return width_; }

int Particle::FindParticle(const char *particleName) {
  for (int i{0}; i != numParticleType_; ++i) {
    if (particleName == particleTypeArray_[i]->GetName()) {
      return i;
    }
  }
  std::cout << "Particella non trovata";
  return -1;
}

const int Particle::GetIndex() const { return index_; };

void Particle::AddParticleType(const char *name, const double mass,
                               const int charge, double width) {
  // contro
  if (numParticleType_ < maxNumParticleType_) {
    if (width = 0.) {
      // auto particleType = new ParticleType{name, mass, charge};
      // particleTypeArray_[numParticleType_] = particleType;
      particleTypeArray_[numParticleType_] =
          new ParticleType{name, mass, charge};
    } else {
      particleTypeArray_[numParticleType_] =
          new ResonanceType{name, mass, charge, width};
    }
    ++numParticleType_;
  } else
    throw std::runtime_error(
        "Exceeded max number of particle types. Please increase "
        "Particle::maxNumParticleType_");
}

Particle::Particle(char *name, double pX, double pY, double pZ)
    : index_{FindParticle(name)}, pX_{pX}, pY_{pY}, pZ_{pZ} {}

void Particle::SetIndex(const char *particleName) {
  int index = Particle::FindParticle(particleName);
  index_ = index;
}

void Particle::SetIndex(int index) {
  if (index <= numParticleType_) {
    index_ = index;
  } else {
    throw std::runtime_error("Index out of range.");
  }
}

void Particle::PrintParticleTypeArray() {
  for (int i{0}; i != numParticleType_; ++i) {
    std::cout << '[' << i << "] = ";
    particleTypeArray_[i]->Print();
  }
}

void Particle::PrintParticle() const {
  std::cout << "Indice del tipo di particella:" << index_ << '\n';
  std::cout << "Nome particella:" << particleTypeArray_[index_]->GetName()
            << '\n';
  std::cout << "Componente impulso X:" << pX_ << '\n';
  std::cout << "Componente impulso Y:" << pY_ << '\n';
  std::cout << "Componente impulso Z:" << pZ_ << '\n';
}

double Particle::GetMass() const {
  return particleTypeArray_[index_]->GetMass();
}

double Particle::EnergyTot() const {
  double m{GetMass()};
  return std::sqrt(m * m + pX_ * pX_ + pY_ * pY_ + pZ_ * pZ_);
}

double Particle::InvMass(Particle &p2) {
  double E1 = EnergyTot();
  double E2 = p2.EnergyTot();
  double sumPsquare = std::pow(pX_ + p2.pX_, 2) + std::pow(pY_ + p2.pY_, 2) +
                      std::pow(pZ_ + p2.pZ_, 2);
  return std::sqrt((E1 + E2) * (E1 + E2) - sumPsquare);
}

void Particle::SetP(double pX, double pY, double pZ) {
  pX_ = pX;
  pY_ = pY;
  pZ_ = pZ;
}
