
# change paramters in here

numParticles = 1000;       # how many particles to use                                                 
times = 10000;              # how many times for one experiment
runs=100;                   # how many experiments runs
maxX=500;                 # width
maxY=260;                 # height
side=1;                   # 1 means stick to one side; 2 means stick to two sides
intervals=20;             # interval change for each series of experiment
upperbounce=280;          # highest height

# Function for particle movement

moveP3 = function(p, dt = .2, directionChange = .2,  
maxX = maxX, maxY = maxY){
           theta = runif(1, 0, 2*pi);
           pout = p;
           
           #20% possibility to change direction 
           if(runif(1,0,1) < directionChange){
              s = sqrt(pout[3]^2 + pout[4]^2); 
              pout[3] = cos(theta)*s;
              pout[4] = sin(theta)*s;
           }
           pout[1] = pout[1] + pout[3]*dt;
           pout[2] = pout[2] + pout[4]*dt;

           #hit the left sides
           if(pout[1] < 0 ){
              pout[3] = -pout[3];
              pout[1] = -pout[1];
              }

          #hit the right sides
           if(pout[1] > maxX){
              pout[3] = -pout[3];
              pout[1] = maxX-(pout[1]-maxX);
         } 

          #hit the top if side=1
          
           if(side==1){
            if( pout[2] > maxY  ){
             pout[4] = -pout[4];
              pout[2] = maxY-(pout[2]-maxY);
             }
           }

          #stick to the top if side=2
          if(side==2){
             if( pout[2]>= maxY){
               if(pout[4]!=0){  
                 xWhereHit = p[1] + (maxY-p[2])/pout[4]*pout[3];
                 yWhereHit = p[2] + (maxY-p[2])/pout[4]*pout[4];
                 pout[1] = xWhereHit;
                 pout[2] = maxY;
                 pout[3] = 0;
                 pout[4] = 0;
                 pout[5] = 1; 
                 pout[6] =  xWhereHit ;      
                 pout[7] =  yWhereHit ;}
            }
           }
 
           ##stick to the bottom 
           if( pout[2] <= 0  ){ 
             if(pout[4]!=0){ 
                 xWhereHit = p[1] + (-p[2]/pout[4])*pout[3];
                 yWhereHit = p[2] + (-p[2]/pout[4])*pout[4];
                 pout[1] = xWhereHit;
                 pout[2] = yWhereHit;
                 pout[3] = 0;
                 pout[4] = 0;
                 pout[5] = 1; 
                 pout[6] =  xWhereHit ;      
                 pout[7] =  yWhereHit ;}
            } 
        return(pout);} 

##how many series of experiment
list=0;
while(maxY<=upperbounce){


#how many experiments
k=maxY;  
for(i in 1: runs){

#creat the first particles
   xcord = sample(1: maxX-1, 1);
   ycord = sample(1:maxY-1, 1);
   theta = runif(1, 0, 2*pi);
   xvel =  6*cos(theta);
   yvel =  6*sin(theta);
   st = 0;
   xwh = -9;
   ywh = -9;
   pm = c(xcord, ycord, xvel, yvel, st, xwh, ywh); # x, y, vx, vy stuck xwh ywh

#create the rest particles
 for(i in 1:(numParticles-1)){
  xcord = sample(1: maxX-1, 1);
  ycord = sample(1:maxY-1, 1);
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
 klist=0;
for (i in 1:times){
    pm = t(apply(pm, 1,function(pin){moveP3(p = pin, maxX = maxX,maxY=maxY);}))
    ilist = c(ilist,i);
    clist = c(clist,sum(pm[,5]));
    par(mfrow = c(1,2), pch = 19);
    plot(pm[,1], pm[,2], col = "blue", xlim = c(0,maxX), ylim = c(0,maxY),
            xlab = "", ylab = ""); 
    plot( ilist, clist, 
          xlim = c(0,times), ylim = c(0,numParticles),
          xlab = "time", ylab = "particles removed",  
          col = "red", grid());  
  }

 # estimate first order models
 # First order model y' = k(numParticles - y) solution: y(x) = numParticles(1- exp(-k*x)
 x = ilist; y = clist;
 df <- data.frame(x, y);
 mCondensationFirstOrd <- nls(y ~ (numParticles*(1-exp(-k*x))),
 data = df,start = list(k = .1),trace = F); 
 s = ilist;
 lines(s, predict(mCondensationFirstOrd, list(x = s)), col = "black", lwd = 3 );
 kest = mCondensationFirstOrd$m$getPars() # get k estimate
 text(1, numParticles, pos = 4,paste("k = ", toString(kest)));
 
 k=c(k,kest);  #store K value
}
 list=rbind(list,k);  
 maxY=maxY+intervals;
};

write.table(list, "clipboard", sep="\t", row.names=FALSE);   #output into excel

