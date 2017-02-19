#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

// equation of state of the Lennard-Jones fluid
int main(void)
{
  int i,j,k;
  int NumberOfCyles;
  double PressureSum,PressureCount,dudlCount,Pressure;
  double ChemicalPotentialSum,ChemicalPotentialCount,Dummy1,Dummy2,Dummy3;
  double EnergySquaredSum,Sum_dudl,EnergySum,EnergyCount;
  VECTOR pos;
  FILE *FilePtrMovie,*FilePtr,*thermo;
  printf("**************** Mc_Nvt ***************\n");
  // initialize system
  ReadInputData();
  FilePtr=fopen("results.dat","w");
  FilePtrMovie=fopen("movie.pdb","w");
  thermo=fopen("therm.dat","a");
  NumberOfCyles=0;
  PressureSum=0.0;
  PressureCount=0.0;
  ChemicalPotentialSum=0.0;
  ChemicalPotentialCount=0.0;
  EnergySquaredSum=0.0;
  EnergySum=0.0;
  EnergyCount=0;
  Sum_dudl=0.0;
  dudlCount=0.0;
  // total energy of the system
  EnergySystem();
  printf(" Total Energy Initial Configuration: %lf\n",TotalEnergy);
  printf(" Total Virial Initial Configuration: %lf\n",TotalVirial);
  RunningEnergy=TotalEnergy;
  RunningVirial=TotalVirial;
  Runningdudl=Totaldudl;
  // start MC-cycles
  for(i=1;i<=2;i++)
  {
    // i=1 equilibration
    // i=2 production
    if (i==EQUILIBRATION)
    {
      NumberOfCyles=NumberOfEquilibrationCycles;
      if (NumberOfCyles!=0) printf(" Start Equilibration\n");
    }
    else
    {
      if(NumberOfCyles!=0) printf(" Start Production\n");
      NumberOfCyles=NumberOfProductionCycles;
    }

    NumberOfAttempts=0;
    NumberOfAcceptedMoves=0;

    // intialize the subroutine that adjust the maximum displacement
    Adjust();

    for(j=0;j<NumberOfCyles;j++)
    {
      // attempt to displace a particle
      for(k=0;k<NumberOfDisplacementsPerCycle;k++)
        Mcmove();

      if(i==PRODUCTION)
      {
        // sample averages

        if((j%SamplingFrequency)==0)
        {
          Sample(j,RunningEnergy,RunningVirial,&Pressure,FilePtr);
          PressureSum+=Pressure;
          PressureCount+=1.0;
          // start modification question 4
        EnergySquaredSum+=SQR(RunningEnergy);
        EnergySum+=RunningEnergy;
        EnergyCount+=1;
        Sum_dudl+=Runningdudl;
        dudlCount+=1.0;
        //EnergyParticle(j,RunningEnergy,RunningVirial,&Pressure,FilePtr);
         // Heat capacity. You can use variables EnergySum, EnergyCount,
          // and EnergySquaredSum.
          // end modification
          // calculate the chemical potential
          // do 10 trial chains
          // calculate the average of [exp(-Beta*Energy)]
          // you can use the subroutine EnergyParticle for this
          for(k=0;k<10;k++)
          {
/*    EnergyParticle needs as input a random x, y, and z, then
      the particle number, in this case -1, because the particle
      is added. Then the loop over the particles, starting with
      particle 1. It returns the values dummy1: the energy, and
      dummy2: the virial, which we do not need here.
*/
              pos.x=RandomNumber()*Box;
              pos.y=RandomNumber()*Box;
              pos.z=RandomNumber()*Box;
              EnergyParticle(pos,-1,0,&Dummy1,&Dummy2,&Dummy3);
              ChemicalPotentialSum+=exp(-Beta*Dummy1);
              ChemicalPotentialCount+=1.0;
          }
        }
      }

      if((j%20)==0) WritePdb(FilePtrMovie);
      if((j%(NumberOfCyles/5))==0)
      {
        printf("======>> Done %d out of %d\n",j,NumberOfCyles);
        // write intermediate configuration to file
        Store();

        // adjust maximum displacements
        Adjust();
      }
    }
    if(NumberOfCyles!=0)
    {
      if (NumberOfAttempts!=0)
      {
        printf("Number Of Att. To Displ. A Part.  : %d\n",NumberOfAttempts);
        printf("Success: %d (%lf)\n",NumberOfAcceptedMoves,100.0*NumberOfAcceptedMoves/NumberOfAttempts);
      }
      // test total energy
      EnergySystem();
      if(fabs(TotalEnergy-RunningEnergy)>1.0e-6)
        printf("######### Problems Energy ################\n");
      if(fabs(TotalVirial-RunningVirial)>1.0e-6)
        printf("######### Problems Virial ################\n");
      printf(" Total Energy End Of Simulation    : %lf\n",TotalEnergy);
      printf("       Running Energy              : %lf\n",RunningEnergy);
      printf("       Difference                  : %lf\n",TotalEnergy-RunningEnergy);
      printf(" Total Virial End Of Simulation    : %lf\n",TotalVirial);
      printf("       Running Virial              : %lf\n",RunningVirial);
      printf("       Difference                  : %lf\n",TotalVirial-RunningVirial);
      // print chemical potential and pressure
      if(i==PRODUCTION)
      {
        printf("Average Pressure                  : %lf\n",
          PressureSum/PressureCount);
        printf("dudl                   : %lf\n",
          Sum_dudl/dudlCount);
        printf("Chemical Potential                : %lf\n",
           -(1.0/Beta)*
          (log((ChemicalPotentialSum/ChemicalPotentialCount)/(NumberOfParticles/CUBE(Box))))
           -(1.0/Beta)*(2.0*(8.0/3.0)*M_PI*(NumberOfParticles/CUBE(Box))*Epsilon*
           ((1.0/3.0)*pow(Sigma/CutOff,9.0)-pow(Sigma/CutOff,3.0))));
        printf("Heat capacity                     : %lf\n",
          ((EnergySquaredSum/EnergyCount)-SQR(EnergySum/EnergyCount))/SQR(Temp)/NumberOfParticles);

        fprintf(thermo,"%lf\n",
          Sum_dudl/dudlCount);
      }
    }
  }
  Store();
  fclose(FilePtrMovie);
  fclose(FilePtr);
  fclose(thermo);
  return 0;
}
