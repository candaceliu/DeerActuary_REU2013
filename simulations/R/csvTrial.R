# Demonstration on how to write to a csv file on the fly.
#
# Approximates a stochastic integral and writes the results to a file.
#
# Each time it makes an approximation it writes the last data point
# from an euler approximation and from a Milstein approximation.
#

# Set the parameters for the approximation
T = 2.0;                  # The end time
N = 5000;                 # The number of steps to take
alpha = 1.0;              # The value of the noise level
numberTrials = 200;       # The number of trials to run.
filename = "example.csv"; # Set the file name

# Create the file and write out a header
cat(file=filename,"euler,milstein\n",append=FALSE);

# Go through and make a number of  approximations to the integral
# Write the value of each approximation at time T to the file.
for (trial in 1:numberTrials)
  {
    cat("Trial number ",trial,"/",numberTrials,"\n");
    # Set the initial conditions
    mEuler    = 0.0;
    mMilstein = 0.0;
    W         = 0.0;

    dt        = T/N;
    sdt       = sqrt(dt);
    t         = 0.0;

    # plot(0,0,xlim=c(0,2),ylim=c(-5,5));
    # Perform the Milstein approximation.
    # Include extra steps to provide a better approximation of W
    lupe = 1;
    while(lupe <= N)
      {
        # First update the approximation for dw.
        dw = sum(rnorm(1))*sdt;

        # Update the euler approximation
        mEuler = mEuler +
          sqrt(2*alpha)*exp(alpha*t)*cos(sqrt(2*alpha)*W)*dw;
        #points(t+dt,mEuler,pch=1,col=2);

        # Update the Milstein approximation
        mMilstein = mMilstein +
          sqrt(2*alpha)*exp(alpha*t)*cos(sqrt(2*alpha)*W)*dw -
            alpha*exp(alpha*t)*sin(sqrt(2*alpha)*W)*(dw*dw-dt);
        #points(t+dt,mMilstein,pch=2,col=3);

        # Update the values of the random walk and time
        W = W + dw;
        t = lupe*dt;
        #points(t+dt,exp(alpha*t)*sin(sqrt(2*alpha)*W),pch=3,col=1);

        # Increment the loop counter
        lupe = lupe + 1;
      }

    # Append the data to the end of the example
    cat(file=filename,mEuler,mMilstein,append=TRUE,sep=",",fill=TRUE);

  }
