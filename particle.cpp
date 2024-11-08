#include "particle.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

const char *ParticleType::GetName() const { return name_; }

double ParticleType::GetMass() const { return mass_; }

int ParticleType::GetCharge() const { return charge_; }

double ParticleType::GetWidth() const { return 0.; }

void ParticleType::Print() const {
  std::cout << "Name: " << name_ << ", Mass: " << mass_
            << ", Charge: " << charge_ << '\n';
}
void ResonanceType::Print() const {
  ParticleType::Print();
  std::cout << ", Width: " << width_ << '\n';
}

double ResonanceType::GetWidth() const { return width_; }

int Particle::numParticleType_{0};
const ParticleType *Particle::particleTypeArray_[maxNumParticleType_];

int Particle::FindParticle(const char *particleName) {
  for (int i{0}; i < numParticleType_; ++i) {
    if (particleName == particleTypeArray_[i]->GetName()) {
      return i;
    }
  }
  std::cout << "Particella non trovata";
  return -1;
}

int Particle::GetIndex() const { return index_; };

void Particle::AddParticleType(const char *name, const double mass,
                               const int charge, double width) {
  // contro
  if (numParticleType_ < maxNumParticleType_) {
    if (width == 0.) {
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

Particle::Particle() : index_{-1}, pX_{0}, pY_{0}, pZ_{0} {}
Particle::Particle(const char *name, double pX, double pY, double pZ)
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
double Particle::GetCharge() const {
  return particleTypeArray_[index_]->GetCharge();
}

double Particle::EnergyTot() const {
  double m{GetMass()};
  return std::sqrt(m * m + pX_ * pX_ + pY_ * pY_ + pZ_ * pZ_);
}

double Particle::InvMass(const Particle &p2) const {
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

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (GetMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (index_ > -1) {  // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;

    massMot += particleTypeArray_[index_]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf(
        "Decayment cannot be preformed because mass is too low in this "
        "channel\n");
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
            pout * cos(theta));
  dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi),
            -pout * cos(theta));

  double energy = sqrt(pX_ * pX_ + pY_ * pY_ + pZ_ * pZ_ + massMot * massMot);

  double bx = pX_ / energy;
  double by = pY_ / energy;
  double bz = pZ_ / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz) {
  double energy = EnergyTot();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * pX_ + by * pY_ + bz * pZ_;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  pX_ += gamma2 * bp * bx + gamma * bx * energy;
  pY_ += gamma2 * bp * by + gamma * by * energy;
  pZ_ += gamma2 * bp * bz + gamma * bz * energy;
}
