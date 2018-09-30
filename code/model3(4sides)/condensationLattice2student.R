# condensation lattice 3.R
# At the specified heights will do the specified number of runs.
# Then it will plot the values of k, their mean, and error bars. 

runs = 5;
heights = seq(40, 60, by = 10); heights;
hlen = length(heights);

resultk = c(0);
resultkt = c(0);
resultkMatrix =  matrix(, nrow = runs, ncol = length(heights));

####################  Start two big for loops

for(h in 1:length(heights)){
for(r in 1:runs)
{
numParticles = 1000;              # how many particles in box
its = 12000;                      
xmax = heights[h];  
ymax = xmax;  # must be even

# ------ make main title for plots
mainTitle = paste(numParticles, 
                 " particles in a ", 
                 xmax, " by ", 
                 ymax, " units box. Run ", (h-1)*runs + r, " of ", runs*length(heights), " simulations.", sep ="");
mainTitle ;

# -------------------------------------------------------------
# Function Definitions
# ----------------------------------------------------------------
e1234 = rbind(c(1,0), c(-1,0),c(0,1),c(0,-1));   # four directions
e1234notStuck = cbind(e1234, c(0,0,0,0));        # add stuck conditional
e1234notStuck;
# ----------------------------------------------------------------
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
# ----------------------------------------------------------
# populate box with particles randomly choosen
 p1 = c(xmax/2, ymax/2, 0) # x, y, stuck
 pm = p1;
 for(i in 1:numParticles){
   xcord = sample(1:xmax, 1);
   ycord = sample(1:ymax, 1);
   stuck = 0;
   pm = rbind(pm,c(xcord, ycord, stuck));  
 }
# -----------------------------------------------------------
 clist = 0;   # particle counts
 ilist = 0;   # time  
# ---------------------------------------------------------
# main part of simulation
# ----------------------------------------------------------
for (i in 1:its){
    pm = t(apply(pm, 1,function(pin){movePLatice(p = pin,v4x2 = e1234notStuck);}))
    ilist = c(ilist,i);
    clist = c(clist,sum(pm[,3]));
 #   par(mfrow = c(1,2), pch = 19);
 #   plot(pm[,1], pm[,2], col = "blue", xlim = c(0,xmax), ylim = c(0,ymax),
 #           xlab = "", ylab = ""); 
 #   plot( ilist, clist, 
 #         xlim = c(0,its), ylim = c(0,numParticles),
 #         xlab = "time", ylab = "particles removed",  
 #         col = "red", grid());  
 #
  }; # end main part of simulation
#---------------------------------------------
# calculation k theory
#---------------------------------------------
twoymax = 2*ymax;
kt = ((2*sum(seq(1,twoymax-1)*seq(twoymax-1,1)))/(twoymax-1))^(-1);kt;
#--------------------------------------------------
# plot data, regression curve , and theory curve
#--------------------------------------------------
# first plot data
#
    plot( ilist, clist, 
          xlim = c(0,its), ylim = c(0,numParticles),
          xlab = "time", ylab = "particles removed",  
          col = "red", font.main = 1,
          main = mainTitle,
          grid()
    );  
#---------------------------------------------------
# then plot theory curve using kt
#
s = seq(1,its);
lines(s, numParticles*(1-exp(-s*kt)), lwd = 6, col = "cyan");
#-------------------------------------------------------
# do regression calculations
#
# First order model y' = k(numParticles - y) solution: y(x) = numParticles(1- exp(-k*x)
x = ilist; y = clist;
df <- data.frame(x, y);
mCondensationFirstOrd <- nls(y ~ (numParticles*(1-exp(-k*x))), #algorithm = "port",
 data = df,
 start = list(k = .01),
 trace = F); # was T
 #s <- seq(from = 0, to = 300, length = 50);
 s = ilist;
 lines(s, predict(mCondensationFirstOrd, list(x = s)), col = "blue", lwd = 3 );
#----------------------------------------------------------
# summary(mCondensationFirstOrd)  # uncomment if you want details
#---------------------------------------------------------
kreg = mCondensationFirstOrd$m$getPars() # get regression k estimate and then plot it
text(1, .9*numParticles, pos = 4,paste("regression k = ", toString(round(kreg,6)) )); 
# -------------------------------------------------------------
# plot theory k value
text(1, .8*numParticles, pos = 4,paste("theory k = ", toString(round(kt,6)) )); 
#-------- draw legend
legend(its/2,numParticles/2,
  c("data","regression", "theory"),
  lty = c(NA_integer_,1,1),
  lwd = c(3,3, 6),
  pch = c(1, NA_integer_,NA_integer_),
  col=c("red","blue","cyan"),
  bty = "n" );
#--------------------------
# put value of regression k into the vector resultk for storage
resultk[r] = kreg;
#----------------------------------------------
}  # end big for loop "run r"  <------------
#------------------------- 
# store results
resultkt[h] = kt;   # k theory stored 
resultkMatrix[,h] = resultk;   # regression k's (from the runs) stored
resultkMatrix;  
#------------------------- 
}  # end big for loop "h heights" <------------

## ^^^^^ End of the two big for loops ^^^^^^^^^
resultkMatrix
resultkt 
#------------------------ plot k values and confidence intervals
# ----------------------
# find min and maximum k values so we know the range of the plot
yminPlot = min(resultkt, resultkMatrix);yminPlot ;
ymaxPlot = max(resultkt, resultkMatrix);ymaxPlot ;
# ----------------------
resultkt;  # here is the theoretical value of k
# ----------------------
#  make an empty plot (no data) of the correct size
plot(1, type="n", xlab="side of box", ylab="k values", 
        xlim =            c(min(heights), max(heights) ), 
        ylim = c(yminPlot ,ymaxPlot)
         );
# ----------------------
# plot the values of regression k contained in the resultkMatrix
for(h in 1:length(heights)){
  points(heights[h] + 0*resultkMatrix[,h], resultkMatrix[,h]);
  }
# ----------------------
# plot a blue triangle at the values of theory k contained in the resultkt vector
for(h in 1:length(heights)){
  points(heights[h],resultkt[h], col = "blue", pch = 17);
  }
# ----------------------
# plot a red cross at the mean of the regression k's
#for(h in 1:length(heights)){
#   points(heights[h],mean(resultkMatrix[,h]), col = "red", pch = 3, cex = 2);
#   }
# ----------------------
# create a vector avg with the means of the regression k's (duplicates above a little)
for(h in 1:length(heights)){avg[h] = mean(resultkMatrix[,h]);}
avg
# ----------------------
# plot put a solid circle at the means of the regression k's
for(h in 1:length(heights)){
  points(heights, avg, pch = 19);
 }
###----------------------------------------------------------
# calculate confidence intervals 
lowerCI = c(0);
upperCI = c(0);
for(h in 1:length(heights)){
 lowerCI[h] = t.test(resultkMatrix[,h])$conf.int[1];
 upperCI[h] = t.test(resultkMatrix[,h])$conf.int[2];
 }
###----------------------------------------------------------
# plot confidence intervals centerered on means of the regression k's using arrows function
arrows(heights, lowerCI, heights, upperCI, length=0.05, angle=90, code=3)
###----------------------------------------------------------
heights  # heights (heights of boxes - note height = width also)
resultkMatrix  # regression k's each column is for different box height
resultkt  # theoretical k's (each is for different box height)


