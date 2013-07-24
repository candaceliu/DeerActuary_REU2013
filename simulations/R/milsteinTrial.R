# Test of the Milstein Method for a simple DE
#
#  dm = w dw, m(0)=0
#
# True solution is m = 1/2 w^2 - 1/2 t
#
#

# Set the parameters for the approximation
T = 2.0;
N = 5000;
alpha = 1.0;

# Allocate the space
mEuler       = vector(length = N + 1);
mEuler[1]    = 0.0;
mMilstein    = vector(length = N + 1);
mMilstein[1] = 0.0;
W            = vector(length = N + 1);
W[1]         = 0.0;
dt           = T/N;
sdt          = sqrt(dt);
t            = seq(0,T,by=dt);

# Perform the Milstein approximation.
# Include extra steps to provide a better approximation of W
lupe = 1;
while(lupe <= N)
  {
    # First make a better approximation of W.
    dw = sum(rnorm(4))*sdt/2;
    W[lupe+1] = W[lupe] + dw;

    # Update the euler approximation
    mEuler[lupe+1] = mEuler[lupe] +
      sqrt(2*alpha)*exp(alpha*t[lupe])*cos(sqrt(2*alpha)*W[lupe])*dw;

   # Update the Milstein approximation
    mMilstein[lupe+1] = mMilstein[lupe] +
      sqrt(2*alpha)*exp(alpha*t[lupe])*cos(sqrt(2*alpha)*W[lupe])*dw -
        alpha*exp(alpha*t[lupe])*sin(sqrt(2*alpha)*W[lupe])*(dw*dw-dt);

    # Increment the loop counter
    lupe = lupe + 1;
  }

# Plot the true solution and then the approximations
true = exp(alpha*t)*sin(sqrt(2*alpha)*W);
plot(t,true,type='l',col=1);
points(t,mEuler,type="p",pch=1,col=2);
points(t,mMilstein,type="p",pch=2,col=3);
legend(0.0,2.0,c("True","Euler","Milstein"),
       col=c(1,2,3),lty=c(1,0,0),pch=c(-1,1,2))

# Print out the errors in the console.
cat("Error for Euler: ",sqrt(sum((true-mEuler)^2)),"\n",
    "Error for Milstein: ",sqrt(sum((true-mMilstein)^2)),"\n");
text(0,-1,paste("Error for Euler: ",sqrt(sum((true-mEuler)^2))),
     adj=c(-0.1,0));
text(0,-1.4,paste("Error for Milstein: ",sqrt(sum((true-mMilstein)^2))),
     adj=c(-0.1,0));

