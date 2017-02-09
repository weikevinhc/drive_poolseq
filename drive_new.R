##input file. 
table <- cross2.shared #should be .shared.txt
gcov <- cross2.1000 #output of bedtools genomecoverage but binned in 1000kb windows

##dmel-r5 chromosome size
r5.chr <- list("X" = 22422827, "2L" = 23011544, "2R" = 21146708, "3L" = 21146708, "3R" = 27905053, "4" = 1351857)

## reference allele bias at each read depth bin
orig.table <- table
A.table <- table[table$chr %in% c("2L","3R","3L","3R","4"),]
X.table <- table[table$chr %in% c("X"),]
meanX <- mean(X.table[,5]+X.table[,6])
meanA <- mean(A.table[,5]+A.table[,6])
A.split <- unlist(strsplit(as.character(A.table[,4]), split = "/"))
A.split <- A.split[seq(1,length(A.split), by = 2)]
A.ref.logic <- (A.split == A.table[,3])
A.logic <- table$chr %in% c("2L","2R","3L","3R", "4")
X.split <- unlist(strsplit(as.character(X.table[,4]), split = "/"))
X.split <- X.split[seq(1,length(X.split), by = 2)]
X.ref.logic <- (X.split == X.table[,3])


refA.l <- list()
nrefA.l <- list()
refX.l <- list()
nrefX.l <- list()
refA <- A.table[A.ref.logic, 5:6]
nrefA <- A.table[!A.ref.logic, 5:6]
refX <- X.table[X.ref.logic, 5:6]
nrefX <- X.table[!X.ref.logic, 5:6]

for (i in 1:100) {
  a <- refA[refA[,1]+refA[,2] == i,1]
  b <- refA[refA[,1]+refA[,2] == i,2]
  refA.l$mean <- c(refA.l$mean, mean(a/(a+b)))
  refA.l$bot <- c(refA.l$bot, quantile(a/(a+b), probs = c(0.25,0.75))[1])
  refA.l$top <- c(refA.l$top, quantile(a/(a+b), probs = c(0.25,0.75))[2])
  refA.l$amean <- c(refA.l$amean, mean(a))
  refA.l$bmean <- c(refA.l$bmean, mean(b))
  a <- nrefA[nrefA[,1]+nrefA[,2] == i,1]
  b <- nrefA[nrefA[,1]+nrefA[,2] == i,2]
  nrefA.l$mean <- c(nrefA.l$mean, mean(a/(a+b)))
  nrefA.l$bot <- c(nrefA.l$bot, quantile(a/(a+b), probs = c(0.25,0.75))[1])
  nrefA.l$top <- c(nrefA.l$top, quantile(a/(a+b), probs = c(0.25,0.75))[2])
  nrefA.l$amean <- c(nrefA.l$amean, mean(a))
  nrefA.l$bmean <- c(nrefA.l$bmean, mean(b))
  a <- refX[refX[,1]+refX[,2] == i,1]
  b <- refX[refX[,1]+refX[,2] == i,2]
  refX.l$mean <- c(refX.l$mean, mean(a/(a+b)))
  refX.l$bot <- c(refX.l$bot, quantile(a/(a+b), probs = c(0.25,0.75))[1])
  refX.l$top <- c(refX.l$top, quantile(a/(a+b), probs = c(0.25,0.75))[2])
  refX.l$amean <- c(refX.l$amean, mean(a))
  refX.l$bmean <- c(refX.l$bmean, mean(b))
  a <- nrefX[nrefX[,1]+nrefX[,2] == i,1]
  b <- nrefX[nrefX[,1]+nrefX[,2] == i,2]
  nrefX.l$mean <- c(nrefX.l$mean, mean(a/(a+b)))
  nrefX.l$bot <- c(nrefX.l$bot, quantile(a/(a+b), probs = c(0.25,0.75))[1])
  nrefX.l$top <- c(nrefX.l$top, quantile(a/(a+b), probs = c(0.25,0.75))[2])
  nrefX.l$amean <- c(nrefX.l$amean, mean(a))
  nrefX.l$bmean <- c(nrefX.l$bmean, mean(b))
}
c <- refA.l$mean
d <- nrefA.l$mean
nA <- sqrt(c*(1-d)/(d*(1-c)))
rA <- c/(nA*(1-c))
corA <- 1/(1/rA+1)
c <- refX.l$mean
d <- nrefX.l$mean
nX <- sqrt(c*(1-d)/(d*(1-c)))
rX <- c/(nX*(1-c))
range <- 25:75
rangeX <- 20:60

plot(rangeX-0.1, refX.l$mean[rangeX], type = "o", pch = 1, col = "blue", cex = 0.8, cex.axis = 0.6, las = 1,
     xlim= c(rangeX[1], range[length(range)]), ylim = c(0,0.5), xlab = "Read depth", ylab = "P1 freq.")
lines(rangeX+0.1, nrefX.l$mean[rangeX], type = "o", pch = 2, col = "blue", cex = 0.8)
lines(rangeX, (refX.l$amean/nX/(refX.l$amean/nX + refX.l$bmean))[rangeX], type = "o", pch = 20, col = "blue")
lines(range-0.1, refA.l$mean[range], type = "o", pch = 1, cex = 0.8)
lines(range+0.1, nrefA.l$mean[range], type = "o", pch = 2, cex = 0.8)
lines(range, (refA.l$amean/nA/(refA.l$amean/nA + refA.l$bmean))[range], type = "o", pch = 20)
legend("topleft", c("Ref", "Nref", "Corrected"), col = c("black"), bty = "n",
       lty = c(1,1,1), pch = c(1,2,20), cex = 0.8, pt.cex = 1)
abline(v = floor(meanX*2), col = "blue", lty = 2)
abline(v = floor(meanA*2), lty = 2)

## reference bias correction
X.logic <- "X" == table$chr
chrsplit <- unlist(strsplit(as.character(table[X.logic,4]), split = "/"))
chrsplit <- chrsplit[seq(1,length(chrsplit), by = 2)]
ref.logic <- (chrsplit == table[X.logic,3])
dp <- table[X.logic,5] + table[X.logic,6]
dpx2.logic <- table[X.logic,5]+table[X.logic,6] < meanX*2
table[X.logic,5] <- ifelse(dpx2.logic, ifelse(ref.logic, table[X.logic,5], table[X.logic,5] * nX[dp]), table[X.logic,5])
table[X.logic,6] <- ifelse(dpx2.logic, ifelse(ref.logic, table[X.logic,6] * nX[dp], table[X.logic,6]), table[X.logic,6])

for (i in c("2L","2R", "3R","3L","3R","4")) {
  A.logic <- i == table$chr
  chrsplit <- unlist(strsplit(as.character(table[A.logic,4]), split = "/"))
  chrsplit <- chrsplit[seq(1,length(chrsplit), by = 2)]
  ref.logic <- (chrsplit == table[A.logic,3])
  meanA <- mean(table[A.logic,5]+table[A.logic,6])
  dp <- table[A.logic,5] + table[A.logic,6]
  dpx2.logic <- table[A.logic,5]+table[A.logic,6] < meanA*2
  table[A.logic,5] <- ifelse(dpx2.logic, ifelse(ref.logic, table[A.logic,5], table[A.logic,5] * nA[dp]), table[A.logic,5])
  table[A.logic,6] <- ifelse(dpx2.logic, ifelse(ref.logic, table[A.logic,6] * nA[dp], table[A.logic,6]), table[A.logic,6])
}

## filtering out SNPs based on cutoffs

snp.filt <- list()
for (n in names(r5.chr)) {
  a <- subset(table, as.character(table$chr) == n)
  f <- a
  if (n == "X") {
    f <- f[(f[,5] + f[,6] <= meanX*2) & (f[,5] + f[,6] >= 20),]
  } else {
    f <- f[(f[,5] + f[,6] <= meanA*2) & (f[,5] + f[,6] >= 25),]
  }
  snp.filt[[n]] <- f
}

##estimating betabinom rho for X and A
library(VGAM)
A.logic <- table$chr %in% c("2L","3R","3L","3R","4")
X.logic <- table$chr %in% c("X")

ol <- c()
ov <- c()
bv <- c()
pv <- c()
nbrho <- c()
bbv <- c()
bbcoef <- c()
abb <- table[A.logic,5][orig.table[A.logic,5]+orig.table[A.logic,6]>=25]
bbb <- table[A.logic,6][orig.table[A.logic,5]+orig.table[A.logic,6]>=25]
fit <- vglm(cbind(round(abb), round(bbb)) ~ 1, betabinomial, trace = TRUE)
rhoA <- Coef(fit)[2]
abb <- table[A.logic,5][orig.table[A.logic,5]+orig.table[A.logic,6]>=25]
bbb <- table[A.logic,6][orig.table[A.logic,5]+orig.table[A.logic,6]>=25]

for (i in 25:60) {
  a <- table[A.logic,5][orig.table[A.logic,5]+orig.table[A.logic,6]==i]
  b <- table[A.logic,6][orig.table[A.logic,5]+orig.table[A.logic,6]==i]
  ov <- c(ov,var(a))
  
  
  bv <- c(bv,var(rbinom(1e6,i,mean(a/(a+b)))))
  pv <- c(pv,var(rpois(1e6,mean(a))))
}

for (i in 25:60) {
  a <- table[A.logic,5][orig.table[A.logic,5]+orig.table[A.logic,6]==i]
  b <- table[A.logic,6][orig.table[A.logic,5]+orig.table[A.logic,6]==i]
  bbv <- c(bbv, var(rbetabinom(1e6, i, mean(a/(a+b)), bbrho )))
}

plot(hdp$breaks[-1],hdp$counts,xlim = c(25,60),
     xaxt="n",yaxt="n",xlab="",ylab="", , xaxs = "i", ylim = c(0, 12500), type = "n", cex.axis = 0.6, cex.lab = 0.8)
polygon(rep(hdp$breaks, each = 2), c(0,0, rep(hdp$counts,each = 2)[-1],0), col = "gray70", border = NA)
axis(4, cex.axis = 0.8)
abline(v = hdp$breaks, col = "white")
par(new=TRUE)
plot(25:60, ov, type ="o", xaxs = "i", ylim = c(5,30), las = 1, pch = 20,
     cex.axis = 0.8, xlab = "read depth", ylab = "variance",lty=2)
lines(25:60, bv, type ="l", col = "forestgreen", lty = 2)
lines(25:60, pv, type ="l", col = "goldenrod1", lty = 2)
lines(25:60, bbv, type ="l", col = "red", lty = 2)
legend("topright", c("Obs", "Binom", "Pois", "Beta"), col = c("black", "forestgreen", "goldenrod1", "red"), bty = "n",
       lty = c(1,2,2,2), cex = 0.8, pt.cex = c(1,0,0,0), pch = c(20))


##bin and plot allele frequencies from filtered snps.
r5.chr <- list("X" = 22422827, "2L" = 23011544, "2R" = 21146708, "3L" = 21146708, "3R" = 27905053, "4" = 1351857)
window = 200000

abb <- round(snp.filt$X[,5])
bbb <- round(snp.filt$X[,6])
Xfit <- vglm(cbind(abb, bbb) ~ 1, betabinomial, trace = TRUE)
Xrho <- Coef(Xfit)[2]
abb <- round(c(snp.filt$`2L`[,5],snp.filt$`2R`[,5],snp.filt$`3L`[,5],snp.filt$`3R`[,5]))
bbb <- round(c(snp.filt$`2L`[,6],snp.filt$`2R`[,6],snp.filt$`3L`[,6],snp.filt$`3R`[,6]))
Afit <- vglm(cbind(abb, bbb) ~ 1, betabinomial, trace = TRUE)
Arho <- Coef(Afit)[2]

chrend <- c()
A.v <- c()
B.v <- c()
gwwin <- c(0)
snp <- c()
chr_sep <- c()
counter <- 1
chr_mark <- c()

for (n in names(r5.chr)) {
  chr_mark <- c(chr_mark,counter)
  a <- snp.filt[[n]]
  
  win <- c(0)
  win.t <- c(0)
  while (win.t[1] < r5.chr[[n]]) {
    sub.tab <- subset(a, a[,2] >= win.t[1] & a[,2] < (win.t[1] + window))
    win.t <- c(win.t[1] + window, win.t)
    if (nrow(sub.tab) > 0) {
      A.v <- c(A.v, sum(sub.tab[,5]))
      B.v <- c(B.v, sum(sub.tab[,6]))
      snp <- c(snp, nrow(sub.tab))
      win <- c(win.t[1] + window, win)
      counter <- counter + 1
    }
    
  }
  win[1] <- r5.chr[[n]]
  gwwin <- c(win[-length(win)] + gwwin[1], gwwin)
  chrend <- c(gwwin[1], chrend)
  
}
chrend <- c(0,rev(chrend))
gwwin <- rev(gwwin[-1])
binomCI <- sapply(1:(chr_mark[2]-1), 
                    function(x) quantile(replicate(10000, 
                                                   {sum(rbetabinom(snp[x], round((A.v[x]+B.v[x])/snp[x]), A.v[x]/(A.v[x]+B.v[x]), rho = Xrho))}),
                                         prob = c(0.005, 0.995))) ## simulating confidence intervals with betabinom for X windows
binomCI <- c(binomCI, sapply(chr_mark[2]:length(snp), 
                                 function(x) quantile(replicate(10000, 
                                                                {sum(rbetabinom(snp[x], round((A.v[x]+B.v[x])/snp[x]), A.v[x]/(A.v[x]+B.v[x]), rho = Arho))}),
                                                      prob = c(0.005, 0.995)))) ## simulating confidence intervals with betabinom for auto windows
binomlow <- binomCI[rep(c(T,F), length(binomCI)/2)]
binomhigh <- binomCI[rep(c(F,T), length(binomCI)/2)]

plot(gwwin, A.v/(A.v+B.v), pch = 20, ylim = c(0.2,0.4), xaxs = "i", yaxs = "i", cex = 0.4, col = "white",
     cex.axis =0.6, las = 1, xaxt = "n", xlab = "Chr coordinates (Mb)", ylab = "P1 freq.", cex.lab = 0.8)

arrows(gwwin, as.numeric(binomhigh/(A.v+B.v)), gwwin, as.numeric(binomlow/(A.v+B.v)),col=gray(0.8), code = 3, length = 0)
points(gwwin, A.v/(A.v+B.v), pch = 20, cex = 0.4)

abline(v=chrend[-c(1, length(chrend))])
lines(c(0,chrend[2]),c(0.333,0.333), col = "red", )
lines(c(chrend[2], tail(chrend, 1)), c(0.25,0.25), col = "red", )
box(lty = 1)

ticks <- c()
ticks.lab <- c()
for (i in 1:(length(chrend)-1)) {
  ticks <- c(ticks, seq(chrend[i],chrend[i+1], 5000000))
  ticks.lab <- c(ticks.lab, seq(chrend[i],chrend[i+1], 5000000)-chrend[i])
}
axis(side = 1, at= ticks, labels = ticks.lab/1e6, cex.axis = 0.7)



##fitting driver at telomere of 2L based on recombination rate/genetic distance
A.v = c()
B.v = c()
win = c()

chrend = c(0)
for (name in c("2L")) {
  snp.tab <- snp.filt[[name]]
  int = c(0,window)
  snp.v = c()
  while (int[2] < snp.tab[nrow(snp.tab),2]) {
    sub.tab <- subset(snp.tab, snp.tab[,2] > int[1] & snp.tab[,2] < int[2])
    if (nrow(sub.tab) > 0) {
      A.v <- c(A.v, sum(sub.tab[,5]))
      B.v <- c(B.v, sum(sub.tab[,6]))
      snp.v <- c(snp.v, length(sub.tab[,5]))
      win= c(win, int[1] + tail(chrend, 1))
    }
    int <- int + window
    
  }
  chrend <- c(chrend, tail(win,1))
}

x <- win/1e6
recomb <- -0.01*x^3 + 0.20*x^2 + 2.59*x - 1.56 ##genetic location equation from Fiston-Lavier
recomb[x < 0.53] <- 0 ## no recombination at the first 500kb
recomb[which(recomb[recomb > 50][1] == recomb):length(recomb)] <- 50 ##regions with rate of > 50 are deemed unlinked and therefore kept at 50
plot(win/1e06, recomb, pch = 20, xaxs = "i",
     type = "o", cex = 0.6, las = 1, cex.axis = 0.8, ylab = "Genetic Distance (cM)", xlab = "Chromosomal Coordinates (Mb)")

lsq <- c()
af <- A.v/(A.v+B.v)
lswin <- seq(0.25,0.28, by = 0.0001)
for (i in lswin) {
  drf <- ((1-recomb/100)*i*2 + (recomb/100)*(1-i*2))/2
  lsq <- c(lsq, sum((drf - af)^2))
}

i <- lswin[order(lsq)][1]
drf <- ((1-recomb/100)*i*2 + (recomb/100)*(1-i*2))/2
lines(win/1e6, drf, col = "blue", lty = 2, lwd = 2)

plot(lswin, lsq, pch = 20, type = "b", cex = 0.1,
     xlab = "Deviation at origin", ylab = "Sum of squared residuals", las = 1, cex.axis = 0.8)
points(i, sort(lsq)[1], col = "red", pch = 20, cex = 1)
text(i, sort(lsq)[1] + 0.0025, i, col = "red")
plot(win, af, pch = 20, cex = 0.8)
lines(win, drf, col = "blue", lty = 2, lwd = 2)


##fitting driver at centromere of 3 based on recombination rate/genetic distance
A.v = c()
B.v = c()
win = c()

chrend = c(0)
for (name in c("3L", "3R")) {
  snp.tab <- snp.filt[[name]]
  int = c(0,window)
  snp.v = c()
  while (int[2] < snp.tab[nrow(snp.tab),2]) {
    sub.tab <- subset(snp.tab, snp.tab[,2] > int[1] & snp.tab[,2] < int[2])
    if (nrow(sub.tab) > 0) {
      A.v <- c(A.v, sum(sub.tab[,5]))
      B.v <- c(B.v, sum(sub.tab[,6]))
      snp.v <- c(snp.v, length(sub.tab[,5]))
      win= c(win, int[1] + tail(chrend, 1))
    }
    int <- int + window
    
  }
  chrend <- c(chrend, tail(win,1))
}


x <- seq(0, chrend[2], by = 2e5)/1e6
recomb <- -0.006*x^3 + 0.09*x^2 + 2.94*x - 2.90 ## genetic location equation from Fiston-Lavier
recomb[recomb < 0] <- 0
recomb_high <- match(sort(recomb, decreasing = T)[1], recomb)
recomb[recomb_high:length(recomb)] <- recomb[recomb_high] ## as per Fiston-Lavier, the genetic distance equation is prevented from declining

recomb_3 <- recomb[recomb_high]-recomb  ##this inverts the recombination map, since the function assumes loci starting from the left

x <- seq(chrend[2], chrend[3], by = 2e5) - chrend[2]
x <- x/1e6
recomb <- -0.004*x^3 + 0.24*x^2 - 1.63*x + 50.26
recomb_low <- match(sort(recomb)[1], recomb)
recomb_floor <- recomb[recomb_low]
recomb[1:recomb_low] <- recomb_floor
recomb <- recomb - recomb_floor
recomb[recomb > 50] <- 50
recomb_3 <- c(recomb_3, recomb)

plot(win/1e06, recomb_3, pch = 20, xaxs = "i",
     type = "o", cex = 0.6, las = 1, cex.axis = 0.8, ylab = "Genetic Distance (cM)", xaxt = "n", xlab = "Chromosomal Coordinates (Mb)")
axis(side = 1, at= c(seq(0, chrend[2]/1e06, by = 5),seq(chrend[2]/1e06, chrend[3]/1e06, by = 5)), 
     labels = c(seq(0, chrend[2]/1e06, by = 5),seq(0, (chrend[3]-chrend[2])/1e06, by = 5)),
     cex.axis = 0.8)
abline(v=chrend[2]/1e06)


lsq <- c()
af <- A.v/(A.v+B.v)
lswin <- seq(0.25,0.28, by = 0.0001)
for (i in lswin) {
  drf <- ((1-recomb_3/100)*i*2 + (recomb_3/100)*(1-i*2))/2
  lsq <- c(lsq, sum((drf - af)^2))
}

i <- lswin[order(lsq)][1]
drf <- ((1-recomb_3/100)*i*2 + (recomb_3/100)*(1-i*2))/2
lines(c(seq(0, chrend[2], by = 2e5)/1e6, seq(chrend[2], chrend[3], by = 2e5)/1e6), drf, col = "blue", lty = 2, lwd = 2)
plot(lswin, lsq, pch = 20, type = "b", cex = 0.1,
     xlab = "Deviation at origin", ylab = "Sum of squared residuals", las = 1, cex.axis = 0.8)
CI <- sqrt(sort(lsq)[1]/(length(af)-1)*solve(af %*% af))*(1.96) ## 95% confidence interval
points(i, sort(lsq)[1], col = "red", pch = 20, cex = 1)
points(c(i-CI, i+CI), c(sort(lsq)[1], sort(lsq)[1]) , col = "red", type = "l", lty = 2)
text(i, sort(lsq)[1] + 0.005, i, col = "red")


