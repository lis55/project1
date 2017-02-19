#include <math.h>
#include "system.h"

// calculates the energy of particle i with particles j=jb,NumberOfparticles
void EnergyParticle(VECTOR pos,int i,int jb,double *En,double *Vir, double *dudl)
{
  double r2,Virij,Enij,dudlij;
  double r2i,r6i,lambda;
  VECTOR dr;
  int j;


  Enij=Virij=dudlij=0.0;
  for(j=jb;j<NumberOfParticles;j++)
  {
    if (j!=i)
    {
      dr.x=pos.x-Positions[j].x;
      dr.y=pos.y-Positions[j].y;
      dr.z=pos.z-Positions[j].z;

      // apply boundary conditions
      dr.x-=Box*rint(dr.x/Box);
      dr.y-=Box*rint(dr.y/Box);
      dr.z-=Box*rint(dr.z/Box);

      r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      // calculate the energy
      if(r2<SQR(CutOff))
      {
        lambda=0.7;
        r2i=1.0/r2;
        r6i=CUBE(r2i);
        Enij+=4.0*(SQR(r6i)-r6i);

        if ((j==NumberOfParticles-1)||(i==NumberOfParticles-1))
        {
        dudlij+=4.0*((5.0*pow(lambda,4.0))*SQR(r6i)-(3.0*pow(lambda,2.0))*r6i);
        }
        // start modification virial
        Virij+=24.0*(2.0*(pow(lambda,5.0))*SQR(r6i)-(pow(lambda,3.0))*r6i);

        // end modification virial
      }
    }
  }
  *En=Enij;
  *dudl=dudlij;
  *Vir=Virij;
}

// calculates total system energy
void EnergySystem(void)
{
  double Eni,Viri, dudlij;
  int i;
  Totaldudl=0.0;
  TotalEnergy=0.0;
  TotalVirial=0.0;
  for(i=0;i<NumberOfParticles-1;i++)
  {
    EnergyParticle(Positions[i],i,i+1,&Eni,&Viri,&dudlij);
    TotalEnergy+=Eni;
    Totaldudl+=dudlij;
    TotalVirial+=Viri;
  }

  // add tail-correction
 TotalEnergy+=NumberOfParticles*(8.0/3.0)*M_PI*(NumberOfParticles/CUBE(Box))*Epsilon*
        ((1.0/3.0)*pow(Sigma/CutOff,9.0)-pow(Sigma/CutOff,3.0));
}
