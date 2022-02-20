source("../R/F_Functions.R")


marg.pdf <- function(sSq, df, alpha, beta)
  {
    ## alpha, beta is the hyper parameter of the inverse gamma distribution
    temp <- (df/2)^(df/2) * gamma( df/2 + alpha )* beta^alpha /gamma( df/2 )/gamma(alpha)
    temp <- temp * (sSq)^(df/2-1)/( df*sSq/2+beta)^(df/2+alpha)

    temp
  }


Dataset="Leukemia"
## Dataset="Colon"

## Leukemia Data Set
if(Dataset=="Leukemia")
  {
    leukemiadata <- read.table("../realdata/Textdata/Leukemia_x.txt", quote="\"", comment.char="")
    class        <- read.delim("../realdata/Textdata/Leukemia_y.txt", header=FALSE) # 1: AML,   2: ALL
    data     <- t(leukemiadata)
  }

## Colon Cancer data
if( Dataset=="Colon")
  {
    coloncancer <- read.table("../realdata/Textdata/colon_x.txt", quote="\"", comment.char="")
    class      <- read.delim("../realdata/Textdata/colon_y.txt", header=FALSE) # 1: tumor,   2: normal
    data       <- t(coloncancer)
  }

n1  <- length(class[class==1]) # n1=47
n2  <- length(class[class==2]) # n2=25
con <- data[ class==1, ]
exp <- data[ class==2, ]

n   <- n1+n2
p   <- ncol(data)  # number of Variables

c0    <- sqrt(1/n1+1/n2)
G     <- ncol(con)  # Dimension
df    <- n1-1 + n2-1        # Digree of freedom
  
s1    <- apply(con, 2, var)
s2    <- apply(exp, 2, var)
sp    <- ( (n1-1)*s1+(n2-1)*s2 )/(n1+n2-2)

sSq   <- (c0^2)*sp
res    <-squeezeVar(sSq, df)
a0 = res$df.prior/2 # prior parameters
b0 = (res$df.prior/2)*res$var.prior

sSq.sort <- sort(sSq)

## For Leukemia data set
if( Dataset=="Leukemia"){
  a0 <- 1.601209
  b0 <- 0.008044375
}

## For Colon data set
if( Dataset=="Colon"){
  a0 = res$df.prior/2 # prior parameters
  b0 = (res$df.prior/2)*res$var.prior
}

marg.cdf <- array(0, length(sSq) )
for( i in 1:length(sSq) )
  {
    marg.cdf[i] <- integrate( marg.pdf, 0, sSq.sort[i], df, a0, b0)$value
  }
pvalue=ks.test( marg.cdf, 'punif')$p.value


## a0s <- a0 + c(-50:50)/500
## b0s <- b0 + c(-50:50)/50000

## pvalue.all <- array(0, c( length(a0s), length(b0s) ) )

## for(ii in 1:length(a0s))
##   for(jj in 1:length(b0s))
##   {
    
##     marg.cdf <- array(0, length(sSq) )
##     for( i in 1:length(sSq) )
##       {
##         marg.cdf[i] <- integrate( marg.pdf, 0, sSq.sort[i], df, a0s[ii], b0s[jj])$value
##       }
##     pvalue=ks.test( marg.cdf, 'punif')$p.value

##     pvalue.all[ii, jj] <- pvalue
    
##     if(pvalue>0.01)
##       print( paste("a0=",a0,", b0=",b0,", pvalue=",pvalue ) )
##   }
