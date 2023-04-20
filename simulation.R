rm(list=ls(all=TRUE))
oneSample <- function(Niter = 100, SampleSize = 16000, ResistantRate = 0.2,
                      KillRateRes = 0,
                      KillRateSens = 1,
                      DropRateBt = 0,
                      DropRateAt = 0){
  # Set the initial labeled cell sample (singleton and twins)
  # Nsingleton <- (SampleSize - 2 * NtwinPairs)
  Nresistant <- SampleSize * ResistantRate
  
  # Label for resistant cells
  sample.labels <- paste0("R", seq(1, Nresistant))
  
  # Append labels for sensitive cells
  sample.labels <- c(
    sample.labels,
    paste0("S", seq(1, (SampleSize - Nresistant)))
  )
  
  BTAT.labels <- rep(0, Niter)
  unique.label.AT <- rep(0, Niter)
  unique.label.BT <- rep(0, Niter)
  preR.cell <- rep(0, Niter)
  preS.cell <- rep(0, Niter)
  preRtrueR <- rep(0, Niter)
  preRtrueS <- rep(0, Niter)
  preStrueR <- rep(0, Niter)
  preStrueS <- rep(0, Niter)
  preStrueSpro <- rep(0, Niter)
  preStrueRpro <- rep(0, Niter)
  preRtrueSpro <- rep(0, Niter)
  preRtrueRpro <- rep(0, Niter)
  
  trueNegRate <- rep(0, Niter)
  truePosRate <- rep(0, Niter)
  FOR <- rep(0, Niter)
  BA <- rep(0, Niter)
  
  for(i in 1:Niter){
    # cat(i,"\r")
    set.seed(i)
    
    # Double all the labels
    labels <- rep(sample.labels, times = 2)
    cells <- seq(1, (SampleSize * 2))
    
    # split labeled cells into preBT and preAT
    preBT <- sample(cells, SampleSize, replace = F)
    preAT <- setdiff(cells, preBT)
    preBT.label <- labels[preBT]
    preAT.label <- labels[preAT]
    
    # Randomly drop some cells from sample
    BT.label <- sample(preBT.label, SampleSize * (1 - DropRateBt))
    unique.label.BT[i] = length(unique(BT.label))
    
    # AT treated by drugs and KillRate of cells
    AT.res <- preAT.label[startsWith(preAT.label, prefix = "R")]
    AT.res <- sample(AT.res, length(AT.res) * (1 - KillRateRes)*(1-DropRateAt))
    AT.sens <- preAT.label[startsWith(preAT.label, prefix = "S")]
    AT.sens <- sample(AT.sens, length(AT.sens) * (1 - KillRateSens)*(1-DropRateAt))
    AT.label <- c(AT.res, AT.sens)
    
    BTAT.labels[i] <- length(intersect(BT.label, AT.label))
    
    # unique.label.AT[i] <- length(unique(AT.label)) - BTAT.labels[i]
    unique.label.AT[i] <- length(unique(AT.label)) 
    
    # Number of pre-resistant, pre-sensitive, resistant, new cells
    preR.cell[i] <- sum(BT.label %in% intersect(BT.label, AT.label))
    preS.cell[i] <- length(BT.label) - preR.cell[i]
    preRtrueR[i] <- sum(startsWith(intersect(BT.label, AT.label), "R"))
    preRtrueS[i] <- sum(startsWith(intersect(BT.label, AT.label), "S"))
    preRtrueRpro[i] <- preRtrueR[i]/(preRtrueS[i] + preRtrueR[i]) * 100
    preRtrueSpro[i] <- preRtrueS[i]/(preRtrueS[i] + preRtrueR[i]) * 100
    
    preStrueR[i] <- sum(startsWith(setdiff(BT.label, AT.label), "R"))
    preStrueS[i] <- sum(startsWith(setdiff(BT.label, AT.label), "S"))
    
    preStrueRpro[i] <- preStrueR[i]/(preStrueS[i] + preStrueR[i]) * 100
    preStrueSpro[i] <- preStrueS[i]/(preStrueS[i] + preStrueR[i]) * 100
    
    # TNR, TPR, and FOR
    trueNegRate[i] <- preStrueS[i]/(preRtrueS[i] + preStrueS[i]) * 100
    truePosRate[i] <- preRtrueR[i]/(preRtrueR[i] + preStrueR[i]) * 100
    FOR[i] <- preStrueR[i]/(preStrueR[i] + preStrueS[i]) * 100
    BA[i] <- (trueNegRate[i] + truePosRate[i])/2
  }
  return(data.frame(
    preRtrueS = preRtrueS, preRtrueR = preRtrueR,
    preStrueS = preStrueS,preStrueR = preStrueR,
    preStrueRpro = preStrueRpro, preStrueSpro = preStrueSpro,
    preRtrueRpro = preRtrueRpro, preRtrueSpro = preRtrueSpro, 
    labelBTAT = BTAT.labels, labelAT = unique.label.AT, labelBT = unique.label.BT,
    cellPreR = preR.cell, cellPreS = preS.cell, TNR = trueNegRate, TPR = truePosRate, FOR = FOR, BA = BA ))
}

library(reshape2)
par1 = seq(0.01, 0.32, 0.03) # Resistant Rate
par2 = seq(0, 0.9, 0.1) # dropout

N = 16000
niter = 100
res.table = data.frame(i = NA, j = NA, TPR = NA, TPR.sd = NA, FOR = NA, FOR.sd = NA, BA = NA, BA.sd = NA)
for(i in 1:length(par1)){
  cat(i,"\r")
  for(j in 1:length(par2)){
    t <- oneSample(
        Niter = niter, 
        SampleSize = N,
        ResistantRate = par1[i],
        DropRateBt = par2[j])
    res = apply(t, 2, mean)
    # # # olaparib 1
    # # dist = sqrt((res["cellPreR"]-528)^2 + (res["cellPreS"]-4181)^2 + (res["labelAT"]-3111)^2 + (res["labelBT"]-4244)^2 + (res["labelBTAT"]-401)^2)
    # # 
    # # # olaparib 2
    # # dist = sqrt((res["cellPreR"]-376)^2 + (res["cellPreS"]-4309)^2 + (res["labelAT"]-930)^2 + (res["labelBT"]-4168)^2 + (res["labelBTAT"]-190)^2)
    # 
    # # carbo 1
    # dist = sqrt((res["cellPreR"]-360)^2 + (res["cellPreS"]-2094)^2 + (res["labelAT"]-4134)^2 + (res["labelBT"]-2281)^2 + (res["labelBTAT"]-296)^2)
    
    res2 = apply(t, 2, sd)
    
    res.table = rbind(res.table, data.frame(i = par1[i], j = par2[j], TPR = res["TPR"], TPR.sd = res2["TPR"], 
                                            FOR = res["FOR"], FOR.sd = res2["FOR"],
                                            BA = res["BA"], BA.sd = res2["BA"]))
  }
}

res.table
res.table = res.table[-1,]

res.table = round(res.table, 2)
output_dir <- "./data/output/simulation/Jing_code/"

if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# setwd("C:\\Users\\jtang\\Dropbox\\FIMM\\Anna\\Data")
res1 = data.frame(acast(res.table, i~j, value.var = "TPR"), check.names = F)
openxlsx::write.xlsx(res1, paste0(output_dir, "res1.xlsx"), rowNames = T)

res2 = data.frame(acast(res.table, i~j, value.var = "FOR"), check.names = F)
openxlsx::write.xlsx(res2, paste0(output_dir, "res2.xlsx"), rowNames = T)

res3 = data.frame(acast(res.table, i~j, value.var = "BA"), check.names = F)
openxlsx::write.xlsx(res3, paste0(output_dir, "res3.xlsx"), rowNames = T)

res4 = data.frame(acast(res.table, i~j, value.var = "TPR.sd"), check.names = F)
openxlsx::write.xlsx(res4, paste0(output_dir, "res4.xlsx"), rowNames = T)

res5 = data.frame(acast(res.table, i~j, value.var = "FOR.sd"), check.names = F)
openxlsx::write.xlsx(res5, paste0(output_dir, "res5.xlsx"), rowNames = T)

res6 = data.frame(acast(res.table, i~j, value.var = "BA.sd"), check.names = F)
openxlsx::write.xlsx(res6, paste0(output_dir, "res6.xlsx"), rowNames = T)


res.table[which(res.table$dist == min(res.table$dist, na.rm = T)),]
t <- oneSample(
  Niter = 100, 
  SampleSize = 10000,
  ResistantRate = 0.2,
  DropRateBt = 0,
  DropRateAt = 0)
apply(t, 2, mean)


preS = 4181
preR = 528


K = 1056
N = 2*K
n = K
E = K*(1-choose(N-N/K,n)/choose(N,n))
E2 = K*(1- choose(N-2,n)/choose(N,n))
E3 = K*(1-(N-n)*(N-n-1)/(N*(N-1)))

K-(K-E3)*N/K # solve to 528

preR = 1:K # preR cells
preR2 = rep(preR, times = 2) # doubling

AT.index = sample(seq_along(preR2), K, replace = F)
AT = preR2[AT.index]
length(unique(AT)) # number of unique labels
length(which(table(AT)==2)) # number of labels with two copies
length(which(table(AT)==1)) # number of labels with one copy, need to be equal to 528

BT.index = setdiff(seq_along(preR2), AT.index)
BT = preR2[BT.index]
length(unique(BT)) # number of unique labels
length(which(table(BT)==2)) # number of labels with two copies
length(which(table(BT)==1)) # number of labels with one copy, need to be equal to 528



ramanujan <- function(n){
  n*log(n) - n + log(n*(1 + 4*n*(1+2*n)))/6 + log(pi)/2
}

nchoosek <- function(n,k){
  factorial(n)/(factorial(k)*factorial(n-k))
} 

bignchoosek <- function(n,k){
  exp(ramanujan(n) - ramanujan(k) - ramanujan(n-k))
}
