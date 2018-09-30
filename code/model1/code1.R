
# change paramters in here

numParticles = 1000;      # how many particles to use                                                 
times = 500;               # how many times for one time
runs=2;                   # how many experiments runs


# Function for particle movement

moveP3 = function(p, dt = .3, directionChange = .2,  
maxX = 100, maxY = 100, edge = .2, tol = 2){
           theta = runif(1, 0, 2*pi);
           pout = p;
           
           ##20% possibility to change direction 
           if(runif(1,0,1) < directionChange){
              s = sqrt(p[3]^2 + p[4]^2); 
              pout[3] = cos(theta)*s;
              pout[4] = sin(theta)*s;
           }
           pout[1] = p[1] + pout[3]*dt;
           pout[2] = p[2] + pout[4]*dt;

           ##hit the left and right sides
           if( (pout[1] < edge) || (pout[1] > (maxX- edge)) ){
              pout[3] = -pout[3];
              pout[1] = p[1] + pout[3]*dt;
              }
            if( pout[2] > (maxY - edge)  ){
             pout[4] = -pout[4];
              pout[2] = p[2] + pout[4]*dt;
             }
 
           ##stick to the bottom
           if( pout[2] <= 0  ){ 
             if(p[4]!=0){ 
                 timeToHit = dt + (p[2]/p[4]);
                 xWhereHit = p[1] + (-p[2]/p[4])*p[3];
                 yWhereHit = p[2] + (-p[2]/p[4])*p[4];
                 #timeLeftAfterHit = -(p[2]/p[4]);
                 pout[1] = xWhereHit;
                 pout[2] = yWhereHit;
                 pout[3] = 0;
                 pout[4] = 0;
                 pout[5] = 1; 
                 pout[6] =  xWhereHit ;      
                 pout[7] =  yWhereHit ;}
            } 
        return(pout);} 


k=0;
for(n in 1:runs){
#creat particles
 p1 = c(50, 50, 6, 0, 0,-9,-9); # x, y, vx, vy stuck xwh ywh
 pm = p1;
 for(i in 1:numParticles){
  xcord = sample(1: maxX-1, 1);
  ycord = sample(1:99, 1);
  theta = runif(1, 0, 2*pi);
  xvel =  6*cos(theta);
  yvel =  6*sin(theta);
   st = 0;
   xwh = -9;
   ywh = -9;
  pm = rbind(pm,c(xcord, ycord, xvel, yvel, st, xwh, ywh))  
 }

#plot graphs
 clist = 0; 
 ilist = 0;
for (i in 1:times){
    pm = t(apply(pm, 1,function(pin){moveP3(p = pin, maxX = 200);}))
    ilist = c(ilist,i)
    clist = c(clist,sum(pm[,5]))
    par(mfrow = c(1,2), pch = 19);
    plot(pm[,1], pm[,2], col = "blue", xlim = c(0,maxX), ylim = c(0,100),
            xlab = "", ylab = ""); 
    plot( ilist, clist, 
          xlim = c(0,times), ylim = c(0,numParticles),
          xlab = "time", ylab = "particles removed",  
          col = "red", grid());  
  };

 #estimate models
 # First order model y' = k(numParticles - y) solution: y(x) = numParticles(1- exp(-k*x)
 x = ilist; y = clist;
 df <- data.frame(x, y);
 mCondensationFirstOrd <- nls(y ~ (numParticles*(1-exp(-k*x))), #algorithm = "port",
 data = df,
 start = list(k = .1),
 trace = F); 
 s = ilist;
 lines(s, predict(mCondensationFirstOrd, list(x = s)), col = "black", lwd = 3 );
 kest = mCondensationFirstOrd$m$getPars() # get k estimate
 text(1, numParticles, pos = 4,paste("k = ", toString(kest) ));
 k=c(k,kest); 
}
