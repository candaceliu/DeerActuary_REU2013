#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Helper functions to set command line options. */
#include <getopt.h>
#include <string.h>

#define DEFAULT_FILE "trial.csv"
#define DEBUG

/* Check to see if the order is right between min and max. */
void checkOrder(double *theMin,double *theMax)
{
  double tmp;
  if(*theMax < *theMin)
    {
      tmp     = *theMax;
      *theMax = *theMin;
      *theMin = tmp;
    }

}

/* Routine to calculate the step size for a given range and number of iterations. */
double calcDelta(double theMin,double theMax,int number)
{
  return((theMax-theMin)/((double)number));
}




/* Routine to generate two normally distributed random numbers */
inline void randNormal(double nu[]) {
  double tmp = sqrt(-2.0*log(drand48()));
  double trig = 2.0*M_PI*drand48();
  nu[0] = tmp*sin(trig);
  nu[1] = tmp*cos(trig);
}





int main(int argc,char **argv)
{

  /* Define the basic run time variables. */
	 double dt;
	 double sdt;
	 double t;
   double initialTime =  0.0;
   double finalTime   = 10.0;
   double m[2];
	 double dW[2];
   int lupe;
	 int timeLupe;
	 int numberIters     = 1000;
   int numberTimeSteps = 25000;

	 /* define some book keeping variables. */
	 double stochasticIntegral;  // The integral used for the sol. to the pop eqn.
	 double z;                   // The transformed (linear) sol. to the pop eqn.
	 double W;                   // The random walk.

	 /* Define the estimated parameters for the problem. */
	 double r1 = 1.7;     // Deer max reproduction rate
	 double h = .16;      // Harvest rate of the deer
	 double F = 28000;    // Carrying capacity of the deer.
	 double rho = .04;    // Bond fund rate of growth: log(1+rate); 
	 double beta = 9;     // cost due to deer collisions .003*3000 */
	 double g = .05;      // Net target rate of growth of the fund.

	 /* Scaled parameter values. */
	 double rtilde;      // scaled growth rate
	 double ftilde;      // scaled carrying capacity
	 double a;           // exp  exponent for sol. to deep eqn.
	 double g0;          // int. constant for deer pop. solution.



   /* define the parameters ranges*/
   double Pmin     = 2000000.0;
   double alphaMin = 0.0;
   double gammaMin = 0.0;

   double Pmax     = 250000.0;
   double alphaMax = 0.1;
   double gammaMax = 0.1;

   double deltaP;
   double deltaAlpha;
   double deltaGamma;

   int numP     = 10;
   int numAlpha = 10;
   int numGamma = 10;
   int lupeP,lupeAlpha,lupeGamma;

   /* define the parameters */
   double P;
   double alpha;
   double gamma;

	 /* Statistical values */
	 float sumX  = 0.0;
	 float sumX2 = 0.0;
	 float sumM  = 0.0;
	 float sumM2 = 0.0;


   /* Define the output parameters. */
   char outFile[1024];
   strcpy(outFile,DEFAULT_FILE);

  /* File stuff */
  FILE *fp;

  /* Set the step values for the parameters. */
  deltaP     = calcDelta(Pmin,Pmax,numP);
  deltaAlpha = calcDelta(alphaMin,alphaMax,numAlpha);
  deltaGamma = calcDelta(gammaMin,gammaMax,numGamma);

  /* Set the number of iterations used in the main loop. */
  dt  = ((finalTime-initialTime)/((float)numberTimeSteps));
	sdt = sqrt(dt);

#ifdef DEBUG
  printf("Starting iteration. %d iterations.\n",numberTimeSteps);
#endif
  fp = fopen(outFile,"w");
  fprintf(fp,"time,P,alpha,gamma,x,m,sumx,sumx2,summ,summ2,N\n");

	/* Set the seed for the random number generator. */
	srand48(time(NULL));


  for(lupeP=0;lupeP<=numP;++lupeP)
    {
      P = Pmin + deltaP*((double)lupeP);

      for(lupeAlpha=0;lupeAlpha<=numAlpha;++lupeAlpha) 
        {
          alpha = alphaMin + deltaAlpha*((double)lupeAlpha);

          for(lupeGamma=0;lupeGamma<=numGamma;++lupeGamma)
            {
              gamma = gammaMin + deltaGamma*((double)lupeGamma);

#ifdef DEBUG
							/* print a notice */
							printf("%f,%f,%f,%f\n",
										 dt*((float)numberTimeSteps),P,alpha,gamma);
#endif

							/* set the scaled parameters */
							rtilde = r1-h;                     // scaled growth rate
							ftilde = (rtilde/r1)*F;            // scaled carrying capacity
							a      = rtilde-0.5*(alpha*alpha); // exp  exponent for sol. to deep eqn.
							g0     = 0.5*alpha*alpha/a;        // int. constant for deer pop. solution.
							// todo - keep track of 1/a. 
							// keep track of exp(*) or factor it out appropriately?

              /* Start the loop. */
							sumX  = 0.0;
							sumX2 = 0.0;
							sumM  = 0.0;
							sumM2 = 0.0;
              for(lupe=0;lupe<numberIters;++lupe)
                {
									/* set the initial conditions. */
									W    = 0.0;
									m[0] = ftilde;
									m[1] = (P-beta*ftilde)/(g-rho);
									stochasticIntegral = 0.0;

									for(timeLupe=0;timeLupe<numberTimeSteps;++timeLupe)
										{
											/* Set the time step. */
											t = ((double)timeLupe)*dt;

											/* Calc. two normally distributed random numbers */
											if(timeLupe%2==0)
													randNormal(dW); // calc. a new set of random numbers.
											else
												dW[0] = dW[1];    // shift the 2nd number into the first slot.
											dW[0] *= sdt;       // scale the change in W to have the proper variance.

											// Update the integral and then update the population and fund balance.
											stochasticIntegral += exp(a*t+alpha*W)*dW[0] +
												0.5*alpha*exp(a*t+alpha*W)*(dW[0]*dW[0]-dt);
											z = rtilde/a - g0*exp(-a*t-alpha*W) -
												((alpha*rtilde)/a)*exp(-a*t-alpha*W)*stochasticIntegral; 
											m[0] = ftilde/z;
											m[1] += (rho*m[1]+P-beta*m[0])*dt - beta*m[0]*dW[0] 
												- 0.5*alpha*beta*m[0]*(dW[0]*dW[0]-dt);

											W += dW[0];
										}

									// Update the tally used for the statistical ensemble
									sumX  += m[0];
									sumX2 += m[0]*m[0];
									sumM  += m[1]*1.0E-1;
									sumM2 += m[1]*m[1]*1.0E-2;
								}

							//  fprintf(fp,"time,P,alpha,gamma,x,m,sumx,sumx2,summ,summ2,N\n");
							fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d\n",
								dt*((float)numberTimeSteps),
											P,alpha,gamma,m[0],m[1],sumX,sumX2,sumM,sumM2,numberIters);


						}

				}

		}

	fclose(fp);
	return(0);
}

