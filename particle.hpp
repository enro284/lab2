#ifndef PARTICLE_HPP
#define PARTICLE_HPP

class ParticleType
{
public:
    ParticleType(const char *name, double mass, int charge)
        : name_{name}, mass_{mass}, charge_{charge} {}

    const char *GetName() const;
    const double GetMass() const;
    const int GetCharge() const;

    virtual void Print() const;

private:
    const char *name_;
    const double mass_;
    const int charge_;
};

class ResonanceType : public ParticleType
{
    const double width_;

public:
    ResonanceType(const char *name, double mass, int charge, double width)
        : ParticleType(name, mass, charge), width_{width} {}

    double GetWidth();
    void Print() const;
};

class Particle
{
private:
    static const int maxNumParticleType_{7};
    static int numParticleType_;
    static const ParticleType *particleTypeArray_[maxNumParticleType_];
    int index_;

    double pX_, pY_, pZ_;

    int FindParticle(const char *particleName_);

public:
    Particle(char *name, double pX = 0., double pY = 0., double pZ = 0.);
    static void AddParticleType(const char *name, double mass, int charge, double width = 0.);
    static void PrintParticleTypeArray();

    int const GetIndex() const;
    void SetIndex(int index);
    void SetIndex(const char *ParticleName);
    void PrintParticle() const;

    double GetMass() const;
    double EnergyTot() const;
    double InvMass(Particle &p);

    void SetP(double pX, double pY, double pZ);
};

#endif
