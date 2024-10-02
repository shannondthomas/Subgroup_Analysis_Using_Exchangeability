################################################################################
# TITLE: distofPNEXCH_summary.R
#
# PURPOSE: Plot PNEXCH distributions from one-arm continuous outcome case
#
# AUTHOR: Shannon Murphy
# DATE CREATED: FEB 22, 2024
################################################################################

## load in data from 
pnexch_distdata <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/pnexch_dist.csv")
pnexch_distdata_unequaln <- read.csv("C:/Users/mushanno/Desktop/Work/Dissertation1/OneArm/ContinuousOutcome/resultsdata/pnexch_dist_unequaln.csv")


############EQUAL SAMPLE

## Plot results
plothistfunc <- function(df, rownum, xbounds = c(0,1), ybounds = c(0,10000)){
  hist(as.numeric(df[rownum,10:10010]), breaks = seq(0,1,0.05),
       main = paste("N=", df$N[rownum],", Var = ", df$v1[rownum], 
                    ", EffDiff = ", df$es2[rownum] - df$es1[rownum],
                    sep = ''),
       xlim = xbounds, ylim = ybounds,
       xlab = "P(Not Exchangeable)")
}

#plot all var = 1
par(mfrow = c(5,6))
plothistfunc(pnexch_distdata, 15*0 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 5, xbounds = c(0,1), ybounds = c(0,10000))

#plot all var = 10
par(mfrow = c(5,6))
plothistfunc(pnexch_distdata, 15*0 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 10, xbounds = c(0,1), ybounds = c(0,10000))


#plot all var = 100
par(mfrow = c(5,6))
plothistfunc(pnexch_distdata, 15*0 + 11, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 11, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 11, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 11, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 11, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 11, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 12, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 12, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 12, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 12, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 12, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 12, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 13, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 13, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 13, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 13, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 13, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 13, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 14, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 14, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 14, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 14, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 14, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 14, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*0 + 15, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*1 + 15, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*2 + 15, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*3 + 15, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*4 + 15, xbounds = c(0.5,1), ybounds = c(0,10000))
plothistfunc(pnexch_distdata, 15*5 + 15, xbounds = c(0.5,1), ybounds = c(0,10000))







##########UNEQUAL SAMPLE

## Plot results
plothistfunc2 <- function(df, rownum, xbounds = c(0,1), ybounds = c(0,10000)){
  hist(as.numeric(df[rownum,10:10010]), breaks = seq(0,1,0.05),
       main = paste("n1=", df$N[rownum]*(1-df$n2[rownum]),
                    ", n2 = ", df$N[rownum]*df$n2[rownum], 
                    ", EffDiff = ", df$es2[rownum] - df$es1[rownum],
                    sep = ''),
       xlim = xbounds, ylim = ybounds,
       xlab = "P(Not Exchangeable)")
}

#plot all n2mult = 0.1
par(mfrow = c(5,6))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 1, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 2, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 3, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 4, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 5, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 5, xbounds = c(0,1), ybounds = c(0,10000))

#plot all n2mult = 0.2
par(mfrow = c(5,6))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 6, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 7, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 8, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 9, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*0 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*1 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*2 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*3 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*4 + 10, xbounds = c(0,1), ybounds = c(0,10000))
plothistfunc2(pnexch_distdata_unequaln, 15*5 + 10, xbounds = c(0,1), ybounds = c(0,10000))


