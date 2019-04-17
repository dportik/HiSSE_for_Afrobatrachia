library(hisse)
library(parallel)
library(ape)
library(phytools)
library(plyr)

#######################################################################
#Access Sean Harrington's useful scripts
#from "Rate heterogeneity across Squamata, misleading ancestral state reconstruction and the importance of proper null model specification"

source("/Users/hyperolius/Dropbox/HiSSE/get_node_state_probs.R")
get_node_state_probs

source("/Users/hyperolius/Dropbox/HiSSE/plot.hisse.states.PIES.R")
plot.hisse.states.PIES

source("/Users/hyperolius/Dropbox/HiSSE/HiSSE.null4.9rate.R")

#######################################################################
#read in tree
phy<-read.nexus("/Users/hyperolius/Dropbox/HiSSE/Afrobatrachia_Tree.nex")

#######################################################################
#read in the data and make dichromatism data a named vector
raw_data <- read.csv("/Users/hyperolius/Dropbox/HiSSE/Dichromatism_BiSSE.csv")
raw_data
raw_data$Species
raw_data$state

states<- data.frame(raw_data$Species, raw_data$state)
states

#######################################################################
#Make q-rate matrices

TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual <- ParDrop(trans.rates.hisse, c(3,5,8,10))

trans.rates.nodual.equalrates <- trans.rates.nodual
trans.rates.nodual.equalrates[!is.na(trans.rates.nodual.equalrates) & !trans.rates.nodual.equalrates == 0] = 1

trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
trans.rates.bisse.equal <- trans.rates.bisse
trans.rates.bisse.equal[!is.na(trans.rates.bisse.equal)] = 1

trans.rates.nodual.threerates <- trans.rates.nodual
to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] <- 1
to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] <- 2
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] <- 3

trans.state1.matrix <- ParDrop(trans.rates.hisse, c(2, 3, 5, 7, 8, 9, 10, 12))
trans.state0.matrix <- ParDrop(trans.rates.hisse, c(3, 6, 5, 8, 9, 10, 11, 12))

#Check rate matrices
trans.rates.bisse
trans.rates.bisse.equal

trans.rates.nodual
trans.rates.nodual.equalrates
trans.rates.nodual.threerates

trans.state1.matrix
trans.state0.matrix

#######################################################################
#Set prior on root

#For user-specified “root.p”, you should specify the probability for each state. If you are doing a hidden model, 
#there will be four states: 0A, 1A, 0B, 1B. So if you wanted to say the root had to be state 0, you would specify 
#“root.p = c(0.5, 0, 0.5, 0)”.
#Here I set a probability of 0 for dichromatism and 1 for monochromatism on the root, with good evidence.

root_prior = c(0,1)
root_prior_hidden = c(0.5, 0, 0.5, 0)

#######################################################################
#Function to run set of 26 models
hisse.fit <- NA
bisse.fit <- NA

RunModel <- function(model.number){
  
  #BiSSE: λ0 λ1 different, ε0 ε1 different, q01 q10 different "FUll"
  if(model.number==1){
    try(hisse.fit <- hisse(phy, states, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse, 
                           root.type = "user", root.p = root_prior))	
  }
  #BiSSE: λ0 λ1 different, ε0=ε1, q01 q10 different
  if(model.number==2){
    try(hisse.fit <- hisse(phy, states, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse, 
                           root.type = "user", root.p = root_prior))
  }
  #BiSSE: λ0=λ1; ε0=ε1, q01 q10 different
  if(model.number==3){
    try(hisse.fit <- hisse(phy, states, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse, 
                           root.type = "user", root.p = root_prior))
  }
  #BiSSE: λ0 λ1 different, ε0 ε1 different, q01=q10
  if(model.number==4){
    try(hisse.fit <- hisse(phy, states, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse.equal, 
                           root.type = "user", root.p = root_prior))	
  }
  #BiSSE: λ0 λ1 different, ε0=ε1, q01=q10
  if(model.number==5){
    try(hisse.fit <- hisse(phy, states, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse.equal, 
                           root.type = "user", root.p = root_prior))
  }
  #BiSSE: λ0=λ1; ε0=ε1, q01=q10
  if(model.number==6){
    try(hisse.fit <- hisse(phy, states, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse.equal, 
                           root.type = "user", root.p = root_prior))
  }
  
  
  #CID-2: λ0A=λ1A; λ0B=λ1B; ε0A=ε1A; ε0B=ε1B; q’s equal
  if(model.number==7){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual.equalrates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #CID-2: λ0A=λ1A; λ0B=λ1B; ε0A=ε1A; ε0B=ε1B; q’s different
  if(model.number==8){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #CID-2: λ0A=λ1A; λ0B=λ1B; ε’s equal, q’s equal
  if(model.number==9){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual.equalrates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #CID-2: λ0A=λ1A; λ0B=λ1B; ε’s equal, q’s different
  if(model.number==10){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  
  
  #CID-4: all λ different; all ε different; q’s equal
  if(model.number==11){
    try(hisse.fit <- hisse.null4(phy, states, trans.type="equal", 
                                 root.type = "user", root.p = root_prior_hidden))
  }
  #CID-4: all λ different; all ε different; q’s three rate 
  if(model.number==12){
    try(hisse.fit <- hisse.null4(phy, states, trans.type="three.rate", 
                                 root.type = "user", root.p = root_prior_hidden))
  }
  #CID-4: all λ different; ε’s equal, q’s equal
  if(model.number==13){
    try(hisse.fit <- hisse.null4(phy, states, eps.anc=rep(1,8), trans.type="equal", 
                                 root.type = "user", root.p = root_prior_hidden))
  }
  #CID-4: all λ different; ε’s equal; q’s three rate 
  if(model.number==14){
    try(hisse.fit <- hisse.null4(phy, states, eps.anc=rep(1,8), trans.type="three.rate", 
                                 root.type = "user", root.p = root_prior_hidden))
  }
  
  
  #HiSSE: all τ different; all ε different; all q’s different
  if(model.number==15){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: all τ different; all ε different;  q’s three rate
  if(model.number==16){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.threerates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: all τ different; all ε different; q’s equal
  if(model.number==17){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.equalrates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  
  
  #HiSSE: all τ different; ε’s equal, all q’s different
  if(model.number==18){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: all τ different; ε’s equal, q’s three rate
  if(model.number==19){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual.threerates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: all τ different; ε’s equal, q’s equal
  if(model.number==20){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual.equalrates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  
  
  #HiSSE: all τ same; all ε same; all q’s different
  if(model.number==21){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,1,1,1), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: all τ same; all ε same;  q’s three rate
  if(model.number==22){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,1,1,1), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual.threerates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: all τ same; all ε same; q’s equal
  if(model.number==23){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,1,1,1), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual.equalrates, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: transition matrix that only allows hidden effects on state 1
  if(model.number==24){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,0,3), eps.anc=c(1,2,0,3), trans.rate=trans.state1.matrix, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #HiSSE: transition matrix that only allows hidden effects on state 0
  if(model.number==25){
    try(hisse.fit <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,0,3), eps.anc=c(1,2,0,3), trans.rate=trans.state0.matrix, 
                           root.type = "user", root.p = root_prior_hidden))
  }
  #CID4 - 9 rates
  if(model.number==26){
    try(hisse.fit <- hisse.null4.mod.9.rates(phy, states, trans.type ="All.no.dual", 
                           root.type = "user", root.p = root_prior_hidden))
  }
  
  save(phy, hisse.fit, file=paste("DichromatismPrior", model.number, "Rsave", sep="."))
}




#######################################################################
#Execute function in parallel to take advantage of 16 cores on this particular machine.
#You'll want to change that if you don't have that many!

mclapply(1:26, RunModel, mc.cores=16)

##############################################################################################################################################
#######################################################################
#Load results for output writing

model1 <- load("DichromatismPrior.1.Rsave")
hisse.fit1 <- hisse.fit
model2 <- load("DichromatismPrior.2.Rsave")
hisse.fit2 <- hisse.fit
model3 <- load("DichromatismPrior.3.Rsave")
hisse.fit3 <- hisse.fit
model4 <- load("DichromatismPrior.4.Rsave")
hisse.fit4 <- hisse.fit
model5 <- load("DichromatismPrior.5.Rsave")
hisse.fit5 <- hisse.fit
model6 <- load("DichromatismPrior.6.Rsave")
hisse.fit6 <- hisse.fit
model7 <- load("DichromatismPrior.7.Rsave")
hisse.fit7 <- hisse.fit
model8 <- load("DichromatismPrior.8.Rsave")
hisse.fit8 <- hisse.fit
model9 <- load("DichromatismPrior.9.Rsave")
hisse.fit9 <- hisse.fit
model10 <- load("DichromatismPrior.10.Rsave")
hisse.fit10 <- hisse.fit
model11 <- load("DichromatismPrior.11.Rsave")
hisse.fit11 <- hisse.fit
model12 <- load("DichromatismPrior.12.Rsave")
hisse.fit12 <- hisse.fit
model13 <- load("DichromatismPrior.13.Rsave")
hisse.fit13 <- hisse.fit
model14 <- load("DichromatismPrior.14.Rsave")
hisse.fit14 <- hisse.fit
model15 <- load("DichromatismPrior.15.Rsave")
hisse.fit15 <- hisse.fit
model16 <- load("DichromatismPrior.16.Rsave")
hisse.fit16 <- hisse.fit
model17 <- load("DichromatismPrior.17.Rsave")
hisse.fit17 <- hisse.fit
model18 <- load("DichromatismPrior.18.Rsave")
hisse.fit18 <- hisse.fit
model19 <- load("DichromatismPrior.19.Rsave")
hisse.fit19 <- hisse.fit
model20 <- load("DichromatismPrior.20.Rsave")
hisse.fit20 <- hisse.fit
model21 <- load("DichromatismPrior.21.Rsave")
hisse.fit21 <- hisse.fit
model22 <- load("DichromatismPrior.22.Rsave")
hisse.fit22 <- hisse.fit
model23 <- load("DichromatismPrior.23.Rsave")
hisse.fit23 <- hisse.fit
model24 <- load("DichromatismPrior.24.Rsave")
hisse.fit24 <- hisse.fit
model25 <- load("DichromatismPrior.25.Rsave")
hisse.fit25 <- hisse.fit
model26 <- load("DichromatismPrior.26.Rsave")
hisse.fit26 <- hisse.fit

#######################################################################
#Summarize results in a output table 

#Model_Matrix <- matrix(ncol=60)
#colnames(Model_Matrix) <- c("Model","loglik","AIC","HiddenStates","parameters", "p", "p", "p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p",
#                            "p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p","p", "p",
#                            "p", "p","p")
#Model_Matrix <- rbind(Model_Matrix, c("model1", hisse.fit1$loglik, hisse.fit1$AIC, hisse.fit1$hidden.states, hisse.fit1$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model2", hisse.fit2$loglik, hisse.fit2$AIC, hisse.fit2$hidden.states, hisse.fit2$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model3", hisse.fit3$loglik, hisse.fit3$AIC, hisse.fit3$hidden.states, hisse.fit3$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model4", hisse.fit4$loglik, hisse.fit4$AIC, hisse.fit4$hidden.states, hisse.fit4$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model5", hisse.fit5$loglik, hisse.fit5$AIC, hisse.fit5$hidden.states, hisse.fit5$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model6", hisse.fit6$loglik, hisse.fit6$AIC, hisse.fit6$hidden.states, hisse.fit6$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model7", hisse.fit7$loglik, hisse.fit7$AIC, hisse.fit7$hidden.states, hisse.fit7$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model8", hisse.fit8$loglik, hisse.fit8$AIC, hisse.fit8$hidden.states, hisse.fit8$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model9", hisse.fit9$loglik, hisse.fit9$AIC, hisse.fit9$hidden.states, hisse.fit9$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model30", hisse.fit10$loglik, hisse.fit10$AIC, hisse.fit10$hidden.states, hisse.fit10$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model11", hisse.fit11$loglik, hisse.fit11$AIC, "TRUE", hisse.fit11$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model12", hisse.fit12$loglik, hisse.fit12$AIC, "TRUE", hisse.fit12$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model13", hisse.fit13$loglik, hisse.fit13$AIC, "TRUE", hisse.fit13$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model14", hisse.fit14$loglik, hisse.fit14$AIC, "TRUE", hisse.fit14$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model15", hisse.fit15$loglik, hisse.fit15$AIC, hisse.fit15$hidden.states, hisse.fit15$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model16", hisse.fit16$loglik, hisse.fit16$AIC, hisse.fit16$hidden.states, hisse.fit16$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model17", hisse.fit17$loglik, hisse.fit17$AIC, hisse.fit17$hidden.states, hisse.fit17$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model18", hisse.fit18$loglik, hisse.fit18$AIC, hisse.fit18$hidden.states, hisse.fit18$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model19", hisse.fit19$loglik, hisse.fit19$AIC, hisse.fit19$hidden.states, hisse.fit19$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model20", hisse.fit20$loglik, hisse.fit20$AIC, hisse.fit20$hidden.states, hisse.fit20$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model21", hisse.fit21$loglik, hisse.fit21$AIC, hisse.fit21$hidden.states, hisse.fit21$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model22", hisse.fit22$loglik, hisse.fit22$AIC, hisse.fit22$hidden.states, hisse.fit22$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model23", hisse.fit23$loglik, hisse.fit23$AIC, hisse.fit23$hidden.states, hisse.fit23$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model24", hisse.fit24$loglik, hisse.fit24$AIC, hisse.fit24$hidden.states, hisse.fit24$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model25", hisse.fit25$loglik, hisse.fit25$AIC, hisse.fit25$hidden.states, hisse.fit25$solution))
#Model_Matrix <- rbind(Model_Matrix, c("model26", hisse.fit26$loglik, hisse.fit26$AIC, hisse.fit26$hidden.states, hisse.fit26$solution))
#Model_Matrix <- Model_Matrix[-1,]
#Model_Matrix
#
#write.table(format(Model_Matrix, digits = 5, scientific=FALSE), file="PriorModel_Summary_Table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#######################################################################
#Do the AIC weights

# AIC weight
aics <- data.frame(c("model1", "model2", "model3", "model4","model5","model6","model7","model8","model9","model10","model11","model12",
                     "model13","model14","model15","model16","model17","model18","model19","model20","model21","model22","model23","model24","model25","model26"), 
                   c(hisse.fit1$loglik,hisse.fit2$loglik,hisse.fit3$loglik,hisse.fit4$loglik,hisse.fit5$loglik,hisse.fit6$loglik,hisse.fit7$loglik,hisse.fit8$loglik,
                     hisse.fit9$loglik,hisse.fit10$loglik,hisse.fit11$loglik,hisse.fit12$loglik,hisse.fit13$loglik,hisse.fit14$loglik,hisse.fit15$loglik,hisse.fit16$loglik,
                     hisse.fit17$loglik,hisse.fit18$loglik,hisse.fit19$loglik,hisse.fit20$loglik,hisse.fit21$loglik,hisse.fit22$loglik,hisse.fit23$loglik,
                     hisse.fit24$loglik,hisse.fit25$loglik,hisse.fit26$loglik),
                   c(hisse.fit1$AIC,hisse.fit2$AIC,hisse.fit3$AIC,hisse.fit4$AIC,hisse.fit5$AIC,hisse.fit6$AIC,hisse.fit7$AIC,hisse.fit8$AIC,
                     hisse.fit9$AIC,hisse.fit10$AIC,hisse.fit11$AIC,hisse.fit12$AIC,hisse.fit13$AIC,hisse.fit14$AIC,hisse.fit15$AIC,hisse.fit16$AIC,
                     hisse.fit17$AIC,hisse.fit18$AIC,hisse.fit19$AIC,hisse.fit20$AIC,hisse.fit21$AIC,hisse.fit22$AIC,hisse.fit23$AIC,hisse.fit24$AIC,hisse.fit25$AIC,hisse.fit26$AIC), 
                   row.names = NULL)
colnames(aics) <- c("model", "loglik", "AIC")
aics
aics <- aics[order(aics$AIC), ]

for(i in 1:dim(aics)[1]){ 
  aics$delta[i] <- aics$AIC[i] - aics$AIC[1]
} 
aics$W <- (exp(-0.5 * aics$delta) / sum(exp(-0.5 * aics$delta)))
aics
write.table(format(aics, digits = 5, scientific=FALSE), file="PriorModel_AIC_Weights.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#######################################################################


#######################################################################
#perform character reconstructions (using the same root prior with monchromatic = 1, dichromatic = 0)

#BiSSE-like set
pp <- MarginRecon(hisse.fit1$phy, hisse.fit1$data, f=hisse.fit1$f, pars=hisse.fit1$solution, hidden.states=hisse.fit1$hidden.states, aic=hisse.fit1$AIC, 
                  root.type = "user", root.p = root_prior, n.cores=16)
save(pp, file="DichromatismPrior.model.1.recon.Rsave")

pp <- MarginRecon(hisse.fit2$phy, hisse.fit2$data, f=hisse.fit2$f, pars=hisse.fit2$solution, hidden.states=hisse.fit2$hidden.states, aic=hisse.fit2$AIC, 
                  root.type = "user", root.p = root_prior, n.cores=16)
save(pp, file="DichromatismPrior.model.2.recon.Rsave")

pp <- MarginRecon(hisse.fit3$phy, hisse.fit3$data, f=hisse.fit3$f, pars=hisse.fit3$solution, hidden.states=hisse.fit3$hidden.states, aic=hisse.fit3$AIC, 
                  root.type = "user", root.p = root_prior, n.cores=16)
save(pp, file="DichromatismPrior.model.3.recon.Rsave")

pp <- MarginRecon(hisse.fit4$phy, hisse.fit4$data, f=hisse.fit4$f, pars=hisse.fit4$solution, hidden.states=hisse.fit4$hidden.states, aic=hisse.fit4$AIC, 
                  root.type = "user", root.p = root_prior, n.cores=16)
save(pp, file="DichromatismPrior.model.4.recon.Rsave")

pp <- MarginRecon(hisse.fit5$phy, hisse.fit5$data, f=hisse.fit5$f, pars=hisse.fit5$solution, hidden.states=hisse.fit5$hidden.states, aic=hisse.fit5$AIC, 
                  root.type = "user", root.p = root_prior, n.cores=16)
save(pp, file="DichromatismPrior.model.5.recon.Rsave")

pp <- MarginRecon(hisse.fit6$phy, hisse.fit6$data, f=hisse.fit6$f, pars=hisse.fit6$solution, hidden.states=hisse.fit6$hidden.states, aic=hisse.fit6$AIC, 
                  root.type = "user", root.p = root_prior, n.cores=16)
save(pp, file="DichromatismPrior.model.6.recon.Rsave")


#Null-2
pp <- MarginRecon(hisse.fit7$phy, hisse.fit7$data, f=hisse.fit7$f, pars=hisse.fit7$solution, hidden.states=hisse.fit7$hidden.states, aic=hisse.fit7$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.7.recon.Rsave")

pp <- MarginRecon(hisse.fit8$phy, hisse.fit8$data, f=hisse.fit8$f, pars=hisse.fit8$solution, hidden.states=hisse.fit8$hidden.states, aic=hisse.fit8$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.8.recon.Rsave")

pp <- MarginRecon(hisse.fit9$phy, hisse.fit9$data, f=hisse.fit9$f, pars=hisse.fit9$solution, hidden.states=hisse.fit9$hidden.states, aic=hisse.fit9$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.9.recon.Rsave")

pp <- MarginRecon(hisse.fit10$phy, hisse.fit10$data, f=hisse.fit10$f, pars=hisse.fit10$solution, hidden.states=hisse.fit10$hidden.states, aic=hisse.fit10$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.10.recon.Rsave")


#Null-4
pp <- MarginRecon(hisse.fit11$phy, hisse.fit11$data, f=hisse.fit11$f, pars=hisse.fit11$solution, aic=hisse.fit11$AIC, 
                  root.type = "user", root.p = root_prior_hidden, four.state.null=TRUE, n.cores=16)
save(pp, file="DichromatismPrior.model.11.recon.Rsave")

pp <- MarginRecon(hisse.fit12$phy, hisse.fit12$data, f=hisse.fit12$f, pars=hisse.fit12$solution, aic=hisse.fit12$AIC, 
                  root.type = "user", root.p = root_prior_hidden, four.state.null=TRUE, n.cores=16)
save(pp, file="DichromatismPrior.model.12.recon.Rsave")

pp <- MarginRecon(hisse.fit13$phy, hisse.fit13$data, f=hisse.fit13$f, pars=hisse.fit13$solution, aic=hisse.fit13$AIC, 
                  root.type = "user", root.p = root_prior_hidden, four.state.null=TRUE, n.cores=16)
save(pp, file="DichromatismPrior.model.13.recon.Rsave")

pp <- MarginRecon(hisse.fit14$phy, hisse.fit14$data, f=hisse.fit14$f, pars=hisse.fit14$solution, aic=hisse.fit14$AIC, 
                  root.type = "user", root.p = root_prior_hidden, four.state.null=TRUE, n.cores=16)
save(pp, file="DichromatismPrior.model.14.recon.Rsave")


#HiSSE set
pp <- MarginRecon(hisse.fit15$phy, hisse.fit15$data, f=hisse.fit15$f, pars=hisse.fit15$solution, hidden.states=hisse.fit15$hidden.states, aic=hisse.fit15$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.15.recon.Rsave")

pp <- MarginRecon(hisse.fit16$phy, hisse.fit16$data, f=hisse.fit16$f, pars=hisse.fit16$solution, hidden.states=hisse.fit16$hidden.states, aic=hisse.fit16$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.16.recon.Rsave")

pp <- MarginRecon(hisse.fit17$phy, hisse.fit17$data, f=hisse.fit17$f, pars=hisse.fit17$solution, hidden.states=hisse.fit17$hidden.states, aic=hisse.fit17$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.17.recon.Rsave")

pp <- MarginRecon(hisse.fit18$phy, hisse.fit18$data, f=hisse.fit18$f, pars=hisse.fit18$solution, hidden.states=hisse.fit18$hidden.states, aic=hisse.fit18$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.18.recon.Rsave")

pp <- MarginRecon(hisse.fit19$phy, hisse.fit19$data, f=hisse.fit19$f, pars=hisse.fit19$solution, hidden.states=hisse.fit19$hidden.states, aic=hisse.fit19$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.19.recon.Rsave")

pp <- MarginRecon(hisse.fit20$phy, hisse.fit20$data, f=hisse.fit20$f, pars=hisse.fit20$solution, hidden.states=hisse.fit20$hidden.states, aic=hisse.fit20$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.20.recon.Rsave")

pp <- MarginRecon(hisse.fit21$phy, hisse.fit21$data, f=hisse.fit21$f, pars=hisse.fit21$solution, hidden.states=hisse.fit21$hidden.states, aic=hisse.fit21$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.21.recon.Rsave")

pp <- MarginRecon(hisse.fit22$phy, hisse.fit22$data, f=hisse.fit22$f, pars=hisse.fit22$solution, hidden.states=hisse.fit22$hidden.states, aic=hisse.fit22$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.22.recon.Rsave")

pp <- MarginRecon(hisse.fit23$phy, hisse.fit23$data, f=hisse.fit23$f, pars=hisse.fit23$solution, hidden.states=hisse.fit23$hidden.states, aic=hisse.fit23$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.23.recon.Rsave")


pp <- MarginRecon(hisse.fit24$phy, hisse.fit24$data, f=hisse.fit24$f, pars=hisse.fit24$solution, hidden.states=hisse.fit24$hidden.states, aic=hisse.fit24$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.24.recon.Rsave")

pp <- MarginRecon(hisse.fit25$phy, hisse.fit25$data, f=hisse.fit25$f, pars=hisse.fit25$solution, hidden.states=hisse.fit25$hidden.states, aic=hisse.fit25$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.25.recon.Rsave")

pp <- MarginRecon(hisse.fit26$phy, hisse.fit26$data, f=hisse.fit26$f, pars=hisse.fit26$solution, hidden.states=TRUE, aic=hisse.fit26$AIC, 
                  root.type = "user", root.p = root_prior_hidden, n.cores=16)
save(pp, file="DichromatismPrior.model.26.recon.Rsave")

######################################################################

rm(hisse.fit1,hisse.fit2,hisse.fit3,hisse.fit4,hisse.fit5,hisse.fit6,hisse.fit7,hisse.fit8,
   hisse.fit9,hisse.fit10,hisse.fit11,hisse.fit12,hisse.fit13,hisse.fit14,hisse.fit15,hisse.fit16,
   hisse.fit17,hisse.fit18,hisse.fit19,hisse.fit20,hisse.fit21,hisse.fit22,hisse.fit23,
   hisse.fit24,hisse.fit25,hisse.fit26)

#######################################################################
#Do the plots using Sean's "plot.hisse.states.PIES" R script

#1 BiSSE: λ0 λ1 different, ε0 ε1 different, q01 q10 different "Full"
model1 <- load("DichromatismPrior.model.1.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#2 BiSSE: λ0 λ1 different, ε0=ε1, q01 q10 different
model2 <- load("DichromatismPrior.model.2.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#3 BiSSE: λ0=λ1; ε0=ε1, q01 q10 different
model3 <- load("DichromatismPrior.model.3.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#4 BiSSE: λ0 λ1 different, ε0 ε1 different, q01=q10
model4 <- load("DichromatismPrior.model.4.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#5 BiSSE: λ0 λ1 different, ε0=ε1, q01=q10
model5 <- load("DichromatismPrior.model.5.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#6 BiSSE: λ0=λ1; ε0=ε1, q01=q10
model6 <- load("DichromatismPrior.model.6.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#7 CID-2: λ0A=λ1A; λ0B=λ1B; ε0A=ε1A; ε0B=ε1B; q’s equal
model7 <- load("DichromatismPrior.model.7.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#8 CID-2: λ0A=λ1A; λ0B=λ1B; ε0A=ε1A; ε0B=ε1B; q’s different
model8 <- load("DichromatismPrior.model.8.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#9 CID-2: λ0A=λ1A; λ0B=λ1B; ε’s equal, q’s equal
model9 <- load("DichromatismPrior.model.9.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#10 CID-2: λ0A=λ1A; λ0B=λ1B; ε’s equal, q’s different
model10 <- load("DichromatismPrior.model.10.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#11 CID-4: all λ different; all ε different; q’s equal
model11 <- load("DichromatismPrior.model.11.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#12 CID-4: all λ different; all ε different; q’s three rate 
model12 <- load("DichromatismPrior.model.12.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#13 CID-4: all λ different; ε’s equal, q’s equal
model13 <- load("DichromatismPrior.model.13.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#14 CID-4: all λ different; ε’s equal, q’s three rate
model14 <- load("DichromatismPrior.model.14.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#15 HiSSE: all τ different; all ε different; all q’s different
model15 <- load("DichromatismPrior.model.15.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#16 HiSSE: all τ different; all ε different;  q’s three rate
model16 <- load("DichromatismPrior.model.16.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#17 HiSSE: all τ different; all ε different; q’s equal
model17 <- load("DichromatismPrior.model.17.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#18 HiSSE: all τ different; ε’s equal, all q’s different
model18 <- load("DichromatismPrior.model.18.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#19 HiSSE: all τ different; ε’s equal, q’s three rate
model19 <- load("DichromatismPrior.model.19.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#20 HiSSE: all τ different; ε’s equal, q’s equal
model20 <- load("DichromatismPrior.model.20.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#21 HiSSE: all τ same; all ε same; all q’s different
model21 <- load("DichromatismPrior.model.21.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#22 HiSSE: all τ same; all ε same;  q’s three rate
model22 <- load("DichromatismPrior.model.22.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#23 HiSSE: all τ same; all ε same; q’s equal
model23 <- load("DichromatismPrior.model.23.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#24 HiSSE: transition matrix that only allows hidden effects on state 1
model24 <- load("DichromatismPrior.model.24.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#25 HiSSE: transition matrix that only allows hidden effects on state 0
model25 <- load("DichromatismPrior.model.25.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)

#26 CID4 - 9 rates
model26 <- load("DichromatismPrior.model.26.recon.Rsave")
plot.hisse.states.PIES(pp, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=4, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)


#######################################################################
#Do a model averaged plot, incorporating all reconstruction results by weight

hisse.list = list()
models <- system(paste("ls -1 ", "DichromatismPrior.model.", "*.recon.Rsave", sep=""), intern=TRUE)
for(i in 1:length(models)){
  load(models[i])
  hisse.list[[i]] = pp
  rm(pp)
}
hisse.list
plot.hisse.states.PIES(hisse.list, rate.param = "net.div", type = "fan", fsize = 0.4, edge.width.rate=3, edge.width.state=2, tip.pie.size=0.2, node.pie.size=0.3)
plot.hisse.states.PIES(hisse.list, rate.param = "net.div", type = "phylogram", fsize = 0.25, edge.width.rate=2, tip.pie.size=0.1, node.pie.size=0.4)


#######################################################################
#Get the support region for values in the best fit model

no19 <- hisse(phy, states, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.nodual.threerates, 
                          root.type = "user", root.p = root_prior_hidden)
no19


SupportRegion(no19, root.type = "user", root.p = root_prior_hidden)
