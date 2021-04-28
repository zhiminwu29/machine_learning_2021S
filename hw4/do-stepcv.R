library(MASS)

## function to do cv with stepAIC
stepcv = function(ddf,yind,xind,fullform,folds,nstep) {

   ##function to extract sse using stepAIC
   keepf = function(mod,maic) {
      yhat = predict(mod,xpred)
      return(sum((yhat-ypred)^2))
   }

   ##null model
   nullform = as.formula(paste(names(ddf)[yind],"~1"))

   ##loop over folds 
   nf = length(unique(folds)) #number of folds
   ssemat = matrix(0,nstep+1,nf)
   for(i in 1:nf) {
      cat("in stepcv on fold: ",i,"\n")

      ypred = ddf[(folds==i),yind]
      xpred = ddf[(folds==i),xind,drop=FALSE]
      n = nrow(ddf)

      nullmod=lm(nullform,ddf[!(folds==i),])

      fwd = stepAIC(nullmod,scope=fullform,direction="forward",
                k=log(n),trace=0,keep=keepf,steps = nstep)

      ssemat[,i]=as.double(fwd$keep)
   }
   return(ssemat)
}

## set up fold id
getfolds = function(nfold,n,dorand=TRUE) {
   fs = floor(n/nfold) # fold size
   fid = rep(1:nfold,rep(fs,nfold))
   diff = n-length(fid)
   if(diff>0) fid=c(1:diff,fid)
   if(dorand) fid = sample(fid,n)
   return(fid)
}

############################################################
if(0) {cat("### simulate data\n")
#simulate data
set.seed(66) # a good seed!!
nsim=1000
psim=20
xsim = matrix(rnorm(nsim*psim),ncol=psim)
bsim=rep(0,psim);bsim[1:3]=1:3
sigma=10.0
ftrue = xsim %*% matrix(bsim,ncol=1)
ysim = ftrue + sigma*rnorm(nsim)
ddfsim = data.frame(y=ysim,x=xsim)
lmall = lm(y~.,ddfsim)
print(summary(lmall))
}

############################################################
if(0) {cat("### try stepcv on the simulate data\n")
lmall = lm(y~.,ddfsim)
nstep=20
set.seed(99)
fid = getfolds(10,nsim)
fcvsim = stepcv(ddfsim,1,2:21,formula(lmall),fid,nstep)

rmse = sqrt(apply(fcvsim,1,sum)/nsim)
shat = summary(lmall)$sigma
yrg = range(c(sigma,shat,rmse))
plot(1:(nrow(fcvsim)-1),rmse[-1], type="b",col="magenta",xlab="numvar",
                        ylab="rmse",cex.axis=1.5,cex.lab=1.5,ylim=yrg)
abline(h=sigma,col="blue",lwd=2)
abline(h=shat,col="red",lwd=2)
}
