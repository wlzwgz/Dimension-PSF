#####testing whether PSF was different from zero(Supplementary Table 1) 
library(metafor)
data<-read_excel("EmpirialData.xlsx", sheet = "PSF-metafor-variance",col_names = TRUE)
data1<-data[which(data$Family=='composite-composite'),]
data2<-data[which(data$Family=='composite-grass'),]
data3<-data[which(data$Family=='composite-legume'),]
data4<-data[which(data$Family=='grass-grass'),]
data5<-data[which(data$Family=='grass-legume'),]
data6<-data[which(data$Family=='legume-legume'),]
simpleanova1<-rma.mv(PSF,
                     V=VarIs_adj,
                     #mods =~ Family,
                     method = "REML",
                     intercept = FALSE,
                     data=data1 )
simpleanova2<-rma.mv(PSF,
                     V=VarIs_adj,
                     #mods =~ Family  ,
                     
                     method = "REML",
                     intercept = FALSE,
                     data=data2 )
simpleanova3<-rma.mv(PSF,
                     V=VarIs_adj,
                     #mods =~ Family  ,
                     
                     method = "REML",
                     intercept = FALSE,
                     data=data3 )
simpleanova4<-rma.mv(PSF,
                     V=VarIs_adj,
                     #mods =~ Family  ,
                     
                     method = "REML",
                     intercept = FALSE,
                     data=data4 )
simpleanova5<-rma.mv(PSF,
                     V=VarIs_adj,
                     #mods =~ Family  ,
                     
                     method = "REML",
                     intercept = FALSE,
                     data=data5 )
simpleanova6<-rma.mv(PSF,
                     V=VarIs_adj,
                     #mods =~ Family  ,
                     
                     method = "REML",
                     intercept = FALSE,
                     data=data6 )
simpleanova1
simpleanova2
simpleanova3
simpleanova4
simpleanova5
simpleanova6

####correlations between microbial dissimilarity and pairwise PSF --Fig. 2e,f and Supplementary Table 2
psf_dis<-read_excel("EmpirialData.xlsx", sheet = "PSF-Microdissimilarity-modeling",col_names = TRUE)
psf_dis <- psf_dis[,c(7:20)]
myd<-corr.test(psf_dis, use = "complete", method = "pearson", adjust = "none")
myd1<-print(myd, digits=3)

####The best model selction--Supplementary Table 3
library(glmulti)
library(MuMIn)
library(randomForest)
library(ggplot2)
library(A3)
library(readxl)

##model selection
data<-read_excel("EmpirialData.xlsx", sheet = "PSF-Microdissimilarity-modeling",col_names = TRUE)
data <- data[1:81,4:16]
data <- data[,-12]#去掉soilrhizobia
global.model <- lm(Pairfeed_total~Soi_nonrhizo + Root_nonrhizo +  Soilpatho + Rootpatho + Soiloomy + Rootoomy +
                     Soilsapro + Rootsapro  + SoilAMF + RootAMF + #Soilrhizobia + 
                     Rootrhizobia , data = data)

summary(global.model)###R-squared 0.1372

####Random forest
set.seed(12345)
train_data<-sample(nrow(data),0.7*nrow(data))
train<-data[train_data,]
test<-data[-train_data,]
fit.forest<-randomForest(Pairfeed_total~.,train,importance=TRUE,mtry=2)
fit.forest
sig <- a3(Pairfeed_total~Soi_nonrhizo + Root_nonrhizo +  Soilpatho + Rootpatho + Soiloomy + Rootoomy +
            Soilsapro + Rootsapro  + SoilAMF + RootAMF + #Soilrhizobia
            Rootrhizobia +0,data, randomForest, p.acc = 0.01, n.folds = 5)

sig    #####R-squared 0.106

####lm is better so use the glmulti model
psf.model<- glmulti(global.model, level = 1, crit = "aicc",confsetsize=2048)

summary(psf.model)
#average all the models
avemodel=data[, c("Soi_nonrhizo" , "Root_nonrhizo" ,  "Soilpatho" , "Rootpatho" , "Soiloomy" , "Rootoomy",
                  "Soilsapro", "Rootsapro", "SoilAMF", "RootAMF", "Rootrhizobia","Pairfeed_total")]
npar=11
modPar=c("Soi_nonrhizo" , "Root_nonrhizo" ,  "Soilpatho" , "Rootpatho" , "Soiloomy" , "Rootoomy",
         "Soilsapro", "Rootsapro", "SoilAMF", "RootAMF",  "Rootrhizobia","Pairfeed_total")
####build a loop
unit=c(1,0)
parEst=rep(unit,each=2^(npar-1))
for (i in 2:npar){
  unit=c(i,0)
  parEst.tmp=rep(rep(unit,each=2^(npar-i)),2^(i-1))
  parEst=cbind(parEst,parEst.tmp)
}
parMat=cbind(parEst[,1:npar],1)
dimnames(parMat)=list(1:(2^npar),modPar)

allModel=list()
for (i in 1:(dim(parMat)[1]-1)) {
  avemodel.tmp=avemodel[,parMat[i,]!=0]
  allModel[[i]]=glm(Pairfeed_total~.,data=avemodel.tmp)
}

modelC=glm(Pairfeed_total~1,data=avemodel)

lm.ave <- model.avg(allModel,modelC)
summary(lm.ave)
#sum of weights
lm.ave$sw


####t-test for whether CE and RYT were greater than 0 and 1-supplementary Table 4

##use plant richness =2 and CE as an example
data<-read_excel("EmpirialData.xlsx", sheet = "PSF-CE-RYT-fig. 3bc.4.S4.S8",col_names = TRUE)
data2 <- data[which(data[,2]==2),]
test <- t.test(data2$CE, mu = 0)
# 显示结果
test 

#####AVOVA test for the dispersion and richness effect--Supplementary Table 5
###use biomass as an example
data<-read_excel("EmpirialData.xlsx", sheet = "PSF-CE-RYT-fig. 3bc.4.S4.S8",col_names = TRUE)
data$SpecRich<-as.factor(data$SpecRich)
data$Dispersion<-as.factor(CM0$Dispersion)
fit <- aov(Biomass.g.m2 ~ SpecRich+Dispersion+Specrich*Dispersion, data=data)
summary(fit)


#####AVOVA test for the dispersion effect on microbial dissimilarity--Supplemantary Table 6
###use soilpatho as an example
data<-read_excel("EmpirialData.xlsx", sheet = "Microbial dissimilarity-fig.2bd",col_names = TRUE)
data$Dispersion<-as.factor(CM0$Dispersion)
fit <- aov(Soilpatho ~ Dispersion, data=data)
summary(fit)


#####AVOVA test for the interaction between PSF and legume--Supplementary Table 7
###use CE as an example
data<-read_excel("EmpirialData.xlsx", sheet = "PSF-CE-RYT-fig. 3bc.4.S4.S8",col_names = TRUE)
data$Legume<-as.factor(data$Legume)
fit <- lm(CE ~ PSF*Legume, data=data)
summary(fit)


####CopulaAnalyses.R
# Load the needed libraries
library(copula)
library(VC2copula)
library(VineCopula)
library(readxl)

# EMPIRICAL DATA
# Load the data and identify species richness levels
TotalData  <- read_excel("EmpirialData.xlsx",col_names = TRUE)
TotalData<-TotalData[,c(2:5)] 
TotalData6 <- TotalData[which(TotalData[,1]==6),]
TotalData3 <- TotalData[which(TotalData[,1]==3),]
TotalData2 <- TotalData[which(TotalData[,1]==2),]

# Test for independence of the variables
mydata <- data.frame()
mydata[1:72,1]    <- -TotalData2[,2]
mydata[1:72,2]    <-  TotalData2[,4]
mydata[73:120,1]  <- -TotalData3[,2]
mydata[73:120,2]  <-  TotalData3[,4]
mydata[121:168,1] <- -TotalData6[,2]
mydata[121:168,2] <-  TotalData6[,4]
cor(mydata, method = "kendall")

# Select the best fitting copula
var_a <- pobs(mydata)[,1]
var_b <- pobs(mydata)[,2]
selectedCopula <- BiCopSelect(var_a, var_b, familyset = 0:10)
selectedCopula

# Test goodness of fit and tail dependency (for one-parameter copula):
gfc <- BiCopGofTest(var_a, var_b,family = selectedCopula$family, par = selectedCopula$par, method = "white", B = 100)
gfc
BiCopPar2TailDep(family = selectedCopula$family, par = selectedCopula$par)

# Test goodness of fit and tail dependency (for two-parameter copula):
gfc <- BiCopGofTest(var_a, var_b,family = selectedCopula$family, par = selectedCopula$par, par2 = selectedCopula$par2, method = "kendall", B = 100)
gfc
BiCopPar2TailDep(family = selectedCopula$family, par = selectedCopula$par, par2 = selectedCopula$par2)

# Perform non-parametric bootstrapping procedure to test tail dependence:
# Create bootstrap replicates (using normal (symmetric) Copula)
dataset<- cbind(var_a,var_b)
BootstrapData <- copsurrog2d(dataset,normalCopula(),corpres = "kendall",numsurrog = 1000)
LB <- seq(0, 0.75, by=0.25)
UB <- seq(0.25, 1, by=0.25)

# Calculate correlation per segment for the actual data:
PartialCorrData <- NA
for (i in 1:length(LB)) {
  datasum <- as.data.frame(dataset[,1]+dataset[,2])
  AboveLower <- which(datasum>2*LB[i])
  BelowUpper <- which(datasum<2*UB[i])
  TC <- intersect(AboveLower,BelowUpper);
  PartialCorrData[i] <- sum((dataset[TC,1]-mean(dataset[,1]))*(dataset[TC,2]-mean(dataset[,2])))/((length(t(datasum))-1)*sqrt(var(dataset[,1])*var(dataset[,2])))
}

# Calculate correlation per segment for the bootstrap replicates:
BSize <-dim(BootstrapData)
PartialCorrBTrap <- data.frame()
for (j in  1:BSize[3]){
  for (i in 1:length(LB)) {
    datasum <- as.data.frame(BootstrapData[,1,j]+BootstrapData[,2,j])
    AboveLower <- which(datasum>(2*LB[i]))
    BelowUpper <- which(datasum<(2*UB[i]))
    TC <- intersect(AboveLower,BelowUpper);
    PartialCorrBTrap[j,i] <- sum((BootstrapData[TC,1,j]-mean(BootstrapData[,1,j]))*(BootstrapData[TC,2,j]-mean(BootstrapData[,2,j])))/((length(t(datasum))-1)*sqrt(var(BootstrapData[,1,j])*var(BootstrapData[,2,j])))
  }
}

# Identify critical threshold for null hypothesis at (1-tailed) 0.05 significance level
Q1 <- sort(PartialCorrBTrap[,1])
Q2 <- sort(PartialCorrBTrap[,2])
Q3 <- sort(PartialCorrBTrap[,3])
Q4 <- sort(PartialCorrBTrap[,4])
LQ1 <- Q1[50]
UQ1 <- Q1[950]
LQ2 <- Q2[50]
UQ2 <- Q2[950]
LQ3 <- Q3[50]
UQ3 <- Q3[950]
LQ4 <- Q4[50]
UQ4 <- Q4[950]

# THEORETICAL DATA
# Load the data
TotalData <- read_excel("AllTheoreticalData.xlsx",col_names = TRUE)

# Fit the best copula for 100 samples with the same structure as the empirical data
FamilySelects <- NA
for (i in 1:100) {
  TotalData6 <- TotalData[which(TotalData[,1]<8),]
  TotalData6 <- TotalData[which(TotalData6[,1]>6),]
  TotalData3 <- TotalData[which(TotalData[,1]<5),]
  TotalData3 <- TotalData[which(TotalData3[,1]>3),]
  TotalData2 <- TotalData[which(TotalData[,1]<4),]
  TotalData2 <- TotalData[which(TotalData2[,1]>2),]
  Datasize6 <- dim(TotalData6)
  Datasize3 <- dim(TotalData3)
  Datasize2 <- dim(TotalData2)
  RandPerm6 <- sample(1:Datasize6[1])
  RandPerm3 <- sample(1:Datasize3[1])
  RandPerm2 <- sample(1:Datasize2[1])
  mydata <- data.frame()
  mydata[1:72,1]     <- -TotalData2[RandPerm2[1:72],2]
  mydata[1:72,2]     <- (TotalData2[RandPerm2[1:72],1]-1)*(TotalData2[RandPerm2[1:72],4]+0.905)
  mydata[73:120,1]   <- -TotalData3[RandPerm3[1:48],2]
  mydata[73:120,2]   <- (TotalData3[RandPerm3[1:48],1]-1)*(TotalData3[RandPerm3[1:48],4]+0.905)
  mydata[121:168,1]  <- -TotalData6[RandPerm6[1:48],2]
  mydata[121:168,2]  <- (TotalData6[RandPerm6[1:48],1]-1)*(TotalData6[RandPerm6[1:48],4]+0.905)
  cor(mydata, method = "kendall")
  var_a <- pobs(mydata)[,1]
  var_b <- pobs(mydata)[,2]
  selectedCopula <- BiCopSelect(var_a, var_b, familyset = NA)
  selectedCopula
  FamilySelects[i] <- selectedCopula$family
}

# Test for independence of the variables for selected example replicate
mydata             <- data.frame()
mydata[1:72,1]     <- -TotalData2[RandPerm2[1:72],2]
mydata[1:72,2]     <- (TotalData2[RandPerm2[1:72],1]-1)*(TotalData2[RandPerm2[1:72],4]+0.905)
mydata[73:120,1]   <- -TotalData3[RandPerm3[1:48],2]
mydata[73:120,2]   <- (TotalData3[RandPerm3[1:48],1]-1)*(TotalData3[RandPerm3[1:48],4]+0.905)
mydata[121:168,1]  <- -TotalData6[RandPerm6[1:48],2]
mydata[121:168,2]  <- (TotalData6[RandPerm6[1:48],1]-1)*(TotalData6[RandPerm6[1:48],4]+0.905)
cor(mydata, method = "kendall")


# Select the best fitting copula
var_a <- pobs(mydata)[,1]
var_b <- pobs(mydata)[,2]
selectedCopula <- BiCopSelect(var_a, var_b, familyset = 0:10)
selectedCopula

# Test goodness of fit and tail dependency (for one-parameter copula):
gfc <- BiCopGofTest(var_a, var_b,family = selectedCopula$family, par = selectedCopula$par, method = "white", B = 100)
gfc
BiCopPar2TailDep(family = selectedCopula$family, par = selectedCopula$par)

# Test goodness of fit and tail dependency (for two-parameter copula):
gfc <- BiCopGofTest(var_a, var_b,family = selectedCopula$family, par = selectedCopula$par, par2 = selectedCopula$par2, method = "kendall", B = 100)
gfc
BiCopPar2TailDep(family = selectedCopula$family, par = selectedCopula$par, par2 = selectedCopula$par2)

# Perform non-parametric bootstrapping procedure to test tail dependence:
# Create bootstrap replicates (using normal (symmetric) Copula)
dataset<- cbind(var_a,var_b)
BootstrapData <- copsurrog2d(dataset,normalCopula(),corpres = "kendall",numsurrog = 1000)
LB <- seq(0, 0.75, by=0.25)
UB <- seq(0.25, 1, by=0.25)

# Calculate correlation per segment for the actual data:
PartialCorrData <- NA
for (i in 1:length(LB)) {
  datasum <- as.data.frame(dataset[,1]+dataset[,2])
  AboveLower <- which(datasum>2*LB[i])
  BelowUpper <- which(datasum<2*UB[i])
  TC <- intersect(AboveLower,BelowUpper);
  PartialCorrData[i] <- sum((dataset[TC,1]-mean(dataset[,1]))*(dataset[TC,2]-mean(dataset[,2])))/((length(t(datasum))-1)*sqrt(var(dataset[,1])*var(dataset[,2])))
}

# Calculate correlation per segment for the bootstrap replicates:
BSize <-dim(BootstrapData)
PartialCorrBTrap <- data.frame()
for (j in  1:BSize[3]){
  for (i in 1:length(LB)) {
    datasum <- as.data.frame(BootstrapData[,1,j]+BootstrapData[,2,j])
    AboveLower <- which(datasum>(2*LB[i]))
    BelowUpper <- which(datasum<(2*UB[i]))
    TC <- intersect(AboveLower,BelowUpper);
    PartialCorrBTrap[j,i] <- sum((BootstrapData[TC,1,j]-mean(BootstrapData[,1,j]))*(BootstrapData[TC,2,j]-mean(BootstrapData[,2,j])))/((length(t(datasum))-1)*sqrt(var(BootstrapData[,1,j])*var(BootstrapData[,2,j])))
  }
}

# Identify critical threshold for null hypothesis at (1-tailed) 0.05 significance level
Q1 <- sort(PartialCorrBTrap[,1])
Q2 <- sort(PartialCorrBTrap[,2])
Q3 <- sort(PartialCorrBTrap[,3])
Q4 <- sort(PartialCorrBTrap[,4])
LQ1 <- Q1[50]
UQ1 <- Q1[950]
LQ2 <- Q2[50]
UQ2 <- Q2[950]
LQ3 <- Q3[50]
UQ3 <- Q3[950]
LQ4 <- Q4[50]
UQ4 <- Q4[950]
