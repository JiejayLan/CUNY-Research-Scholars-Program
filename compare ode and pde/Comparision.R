# runs simulation (no instantanous mixing), 
# compares ODE model to PDE model
# based on: condensation lattice 1.R
#
##############################################################
numParticles = 1000;              # how many particles in box
its = 2000;                      
xmax = 50;  
ymax = 50;  # must be even
#############################################################
# YOU DO NOT NEED TO CHANGE ANYTHING BELOW HERE 
#############################################################
# Function Definitions
#######################################
e1234 = rbind(c(1,0), c(-1,0),c(0,1),c(0,-1));   # four directions
e1234notStuck = cbind(e1234, c(0,0,0,0));        # add stuck conditional
e1234notStuck;
#############################################################
movePLatice = function(p, v4x2, maxX = 100, maxY = 100){
          pout = p;
          if( pout[2] != 0){
           pout = p + v4x2[sample(4,1),];
           if( pout[1] < 0){pout[1] = -pout[1];};
           if( pout[1] > xmax){pout[1] = xmax + (xmax - pout[1]);};
           if( pout[2] > ymax){pout[2] = ymax + (ymax - pout[2]);};
           if(pout[2] ==0){pout[3] = 1};
          }
          return(pout);
          }
####################################### 3
 p1 = c(xmax/2, ymax/2, 0) # x, y, stuck
 pm = p1;
 for(i in 1:numParticles){
   xcord = sample(1:xmax, 1);
   ycord = sample(1:ymax, 1);
   stuck = 0;
   pm = rbind(pm,c(xcord, ycord, stuck));  
 }
####################################### 3
 clist = 0; 
 ilist = 0;
####################################### 3
for (i in 1:its){
    pm = t(apply(pm, 1,function(pin){movePLatice(p = pin,v4x2 = e1234notStuck);}))
    ilist = c(ilist,i);
    clist = c(clist,sum(pm[,3]));
    par(mfrow = c(1,2), pch = 19);
  # ploting slows down simulation so I've commented it out 
  #  plot(pm[,1], pm[,2], col = "blue", xlim = c(0,xmax), ylim = c(0,ymax),
  #          xlab = "", ylab = ""); 
  #  plot( ilist, clist, 
  #        xlim = c(0,its), ylim = c(0,numParticles),
  #        xlab = "time", ylab = "particles removed",  
  #        col = "red", grid());  
  };
 #  plot last "frame" of simulation 
  plot(pm[,1], pm[,2], col = "blue", xlim = c(0,xmax), ylim = c(0,ymax),
          xlab = "", ylab = ""); 
  plot( ilist, clist, 
        xlim = c(0,its), ylim = c(0,numParticles),
        xlab = "time", ylab = "particles removed",  
        col = "red", grid());  


################################# not used
#  twoymax = 2*ymax;
# k from theory 
# kt = ((2*sum(seq(1,twoymax-1)*seq(twoymax-1,1)))/(twoymax-1))^(-1);kt;
# s = seq(1,its);
# lines(s, numParticles*(1-exp(-s*kt)), col = "blue", lwd = 2);
################################# not used

################################################# PDE part
ns = seq(1, 2027, by = 2); 

Ft = function(t){
 k=1/4; 
  L= ymax*2;
  n = ns;
 sum(  ((8*numParticles)/(n^2*pi^2))*exp( -(n^2*pi^2*k*t)/(L^2)  ) );
} 

Ft(0); # should equal numParticles
lines(s,numParticles - sapply(s,Ft), col = "cyan", lwd = 4); #plot PDE sol
#--------------------------------------------------------------------
# do regression calculations for k from ODE
# First order model y' = k(numParticles - y) solution: y(x) = numParticles(1- exp(-k*x)
x = ilist; y = clist;
df <- data.frame(x, y);
mCondensationFirstOrd <- nls(y ~ (numParticles*(1-exp(-k*x))), #algorithm = "port",
 data = df,
 start = list(k = .01),
 trace = F); # was T
 s = ilist;
 lines(s, predict(mCondensationFirstOrd, list(x = s)), col = "black", lwd = 3 );

#############################################################################plot theorical ODE
kt = (1/4)*(ymax)^(-1);kt;
s = seq(1,its);
lines(s, numParticles*(1-exp(-s*kt)), lwd = 6, col = "blue");
#############################################################################plot theorical ODE

#-------- draw legend
legend(.4*its,.45*numParticles,
  c("data (no mixing)","ODE (k from regression)", "PDE (k = 1/4)","ODE from theorical"),
  lty = c(NA_integer_,1,1),
  lwd = c(3,3,6,9),
  pch = c(1, NA_integer_,NA_integer_),
  col=c("red", "black" ,"cyan","blue"),
  bty = "n" );
#--------------------------


