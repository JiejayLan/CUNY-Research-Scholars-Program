#######################################
# Function Def
#######################################
#         Heat/Diffusion open ends
#         x = temp at position x[i]
tsO = function(x){  
 y = x; n = length(x);
 y[1] = 0;  y[n] = 0;
 for(i in 2:(n-1)){
  y[i] = .5*(x[i-1] + x[i+1])
  };
 return(y);}
###################################### open ends 1
particles = 1000;
L = 80;
x = 0:L;
####################################
xone = c();
for(i in 1:L){xone[i] = 1;};
xone;
####################################
uIC = (particles/(L-2))*xone; 
nx = length(x);nx;
x = uIC; 
yM = max(x); ym = 0; 
k = .5;
###################################
its = 6000;
ilist = c();
clist = c();
for (i in 1:its) {
  mainS = paste("time step:",i);
  mainD = paste(particles, "particles");
  ilist = c(ilist,2*i); # double time
    clist = c(clist,sum(x));
    par(mfrow = c(1,2), pch = 19);
      plot( ilist, (particles  - clist), 
          xlim = c(0,2*its), ylim = c(0,particles),
          xlab = "time", ylab = "particles removed",  
          col = "red",
          main = mainD);  
  plot(x, ylim = c(ym,yM), main = mainS, 
          xlab = "x", ylab = "particles");
  x = tsO(x);
};
###################################


