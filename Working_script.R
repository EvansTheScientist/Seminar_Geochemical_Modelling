# This is a combinatorial inverse modeling using groundwater chemical data from the Pra Basin in Ghana
library(RedModRphree)      # R coupled PHREEQC library to run the simulations
library(magrittr)          #provides a set of tools for writing clean and readable code
library(readxl)            #to get data out of Excel and into R
source("rfun.R")           # file that contains all the necessary functions to run the simulations
                           # It has to be in the folder together with other files including, data, and the
                           # working script

### Initial solution coded as a standard PHREEQC input solution
ewat <- c("SOLUTION 1",
          "units mol/kgw",
          "temp 25",
          "pressure 1",
          "pH 6.087",
          "Al 1.038e-07",
          "C(4)  6.802e-04",
          "Ca  1.822e-04",
          "Cl 3.639e-04",
          "Fe   3.349e-07",
          "K  8.952e-05",
          "Mg 1.234e-04",
          "Na 7.743e-04",
          "S(6)  2.041e-05",
          "Si   4.045e-04",
          "PURE 1",
          "END")

## This is the standard "phreeqc.dat" database without unused stuff at
## the end and with added Phlogopite and Plagioclase
phreeqc::phrLoadDatabase("phreeqc_invPra.dat")

## I splitted the groundwater models into 4 single ones: we read them
## scripts, we run them and import the results
fil <- list.files(pattern="Sol.*pqi")
inps <- lapply(fil, RPhreeFile)

phreeqc::phrSetOutputStringsOn(TRUE)
outl <- lapply(inps, function(x) {
  phreeqc::phrRunString(x)
  phreeqc::phrGetOutputStrings()
})

## This code uses the "lapply" function to apply the ReadOut function to each 
# element of the list "outl", and it collects the first element of the result for 
# each element into a new list stored in the variable "res"

res <- lapply(outl, function(x) ReadOut(x)[[1]]) ## what is the meaning of 1 here
#temp <- InputFromList(res[[1]])
## temp contains the input solution in phreeqc script for the starting solution

### Target solution is the second elemental concentration of the res list
target <- res[[1]]$tot
cvec <- target$molal
names(cvec) <- rownames(target)
## We also need pH
cvec <- c(cvec, pH=res[[1]]$desc["pH",])
##cvec <- c(cvec, pH=res[[3]]$desc["pH",])
## Remove the parentheses from names for consistency
names(cvec) <- sub("\\(.\\)","",names(cvec))

## provide a vector of the ion concentrations of the initial solution
init <- c(C=6.802e-04, Ca=1.822e-04, Cl=3.639e-04, Fe=3.349e-07,
          K=8.952e-05, Mg=1.234e-04, Na=7.743e-04, S=2.041e-05, Si=4.045e-04,
          pH=6.087)

## Here we include the primary and secondary mineral phases

primary <- c("Albite", "Plagioclase", "Anorthite", "Phlogopite", "K-mica", "K-feldspar", 
             "Fe(OH)3(a)", "Calcite")
secondary <- c("Kaolinite", "Ca-Montmorillonite", "Chlorite(14A)", "Quartz", "Chalcedony")


tot <- c(DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=3, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=4, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=5, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=6, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=7, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=8, procs=6))

######################  Results
filtered <- Filter(tot, delta=0.5)

#This line of R code utilizes the sapply function to iterate over the names in the 
#object cvec. Within each iteration, it executes the RPinfo function with parameters 
#filtered, ifelse(x=="pH", "desc", "tot"), and x.
#The ifelse statement within the sapply function checks if the current name x is 
#equal to "pH". If true, it uses "desc" as the second argument for RPinfo; otherwise,
#it uses "tot".

dat <- sapply(names(cvec), function(x) RPinfo(filtered, ifelse(x=="pH", "desc", "tot"), x))

## we calculate the rrmse between the modeled and the target concentrations
ind_all <- ComputeMetric(dat, target = cvec, FUN = "rrmse")

ind <- which.min(ind_all) # ind produces the phreeqc calculated list of all parameters

## Inspect the data structure: "filtered" is a list of list 
filtered[[ind]]      ## this just prints it
str(filtered[[ind]]) ## this brings out the structure of it

## Use the function InputFromList() from the RedModRphree package: you
## have documentation!
?InputFromList

## filtered[[ind]] is the best matching combination we found until
## now. From these computed results, form a new PHREEQC input script:
NewInput <- InputFromList(filtered[[ind]])

################################################################################

## Count occurrences of minerals within the top 50 using the results of the rrmse

################################################################################
# ordering the best 50 matched simulations from lowest rrmse to the highest
inds_best <- order(ind_all)[1:50]
# The corresponding rrmse values are:
best50_rrmse <- ind_all[order(ind_all)[1:50]]
## Extract all phase names from the best 50
AllPhases <- unlist(lapply(filtered[inds_best], function(x) rownames(x$pphases)))
#write.csv(AllPhases,"C:\\Users\\Asus\\Desktop\\Marco_combinv-main\\CombInv\\viz\\Frequency.csv", row.names = FALSE)
## Count and visualize the mineral counts
table<-table(AllPhases)
par(mar=c(12.4,5,1,1))
out <- barplot(table, ylab="Frequency",cex.axis=1.6, cex=1.6, cex.lab=1.6, las=2, col=c("cyan"))

##################################################################################

## Visualizing the components of the initial, target, and the simulated
# for the overall best matched simulation

#################################################################################
viz <- rbind(init, dat[which.min(ind_all),], cvec)
# to select the overall best matched simulation we use the "which.min()"
# we then create the barplot of the elemental concentrations of the simulated, initial and the target 
par(mfrow=c(1,1))
out <- barplot(viz, beside=TRUE, log="y", ylim = c(1E-9,10),
               col=c("orange", "light green", "grey"), cex.axis=1.2, cex=1)
mtext(side=2, line=3, "mol/kgw", font=1)

legend("topleft", c("initial","rrmse", "target"),
       fill=c("orange", "light green", "grey"), bty="n", cex=0.8)

#########################################################################

## Plotting the mineral phases of the overall best model

########################################################################
## we apply this code to select the mineral phases contained in the final
## solution of the overall best simulation
best_Sim <- filtered[[which.min(ind_all)]]$pphases
barplot(best_Sim$delta, names.arg=row.names(best_Sim), las=2, cex.axis=1.2, cex=1.3)
## this helps to offset the labels on the axis
mtext(side=2, line=4, "mol/kgw", font=1)

####################################################################################

## To visualize the best 50 matching simulation and the range of the target solution

################################################################################

## Look at argument "na.rm": if it is FALSE, and one of F and K is NA,
## then NA is returned
b1 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", na.rm=FALSE) ##comp=comps,
## dat is a matrix holding the simulation results

## This is the previous behavior (na.rm=TRUE): just compute the RRMSE
## using the only components which are not NA!
b2 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", na.rm=TRUE) ##comp=comps,

viz1 <- dat[order(b1),][1:50,]
viz2 <- dat[order(b2),][1:50,]

## We read in the raw data "Pra_data_M_Good.xlsx" of the samples which the median concentration was taken from.
## You do not need this is you have only one sample representing your initial and final
## composition
samples <- read_excel("Pra_data_M_Good.xlsx", sheet=3)
samples[samples==0] <- NA

## we visualize the simulations and the range of the sample data for the 50 best matched simulations
par(mar=c(2.5,6,1,1))
PlotComb(res=viz1, samples=samples, cex.axis=1.7, cex=1.7, cex.lab=1.4)
mtext(side=2, line=4.5, "mol/kgw", cex=2)

##############################################################################

## Extract the ranges from the "samples" tibble excluding the first 2 columns

#############################################################################
ranges <- sapply(samples[, 3:12], range)
## dat is a matrix holding the simulation results

IsInRange <- function(matsim, ranges) {
  ## We check concentration-wise (column-wise) if we are within the
  ## range. This return a matrix of the same dimension as "matsim"
  ## filled with TRUE, FALSE or NA
  sapply(colnames(ranges), function(conc) ifelse(ranges[1, conc]< matsim[,conc] & ranges[2,conc] > matsim[,conc], TRUE, FALSE))
}

inrange <- IsInRange(matsim=dat, ranges=ranges)

## How to use: TRUE equals 1 and FALSE 0, so we want all the rows
## whose sum is equal the number of columns in "ranges":
num_Inrange1 = which(rowSums(inrange)==ncol(ranges)) ## 0: NO SIMULATION FALLS WITHIN THE RANGE

## Repeat excluding Fe and K (resp. column 4 and 5 in "ranges")
inrange_nokfe <- IsInRange(matsim=dat, ranges=ranges[, -c(4,5)])

num_Inrange2 = which(rowSums(inrange_nokfe)==ncol(ranges[, -c(4,5)])) ## 229

################################################################################

## Here we engage in the process of choosing either individual ions or combinations of 
## ions to assess which simulations best match the target
################################################################################
colnames(dat)
restot <- as.data.frame(lapply(colnames(dat), function(x) ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=x)),USE.NAMES=FALSE)
colnames(restot)<- colnames(dat)

##comps<- c("Ca")
comps<- c("Mg")
#comps<- c("K")
##comps <- c("Ca", "Mg", "Na")
##comps <- c("Na", "K")
#comps <- c("Na", "Ca", "Al")
##comps <- c("Ca", "Mg", "Na")

res_CaMg <- ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=comps)

## Count occurrences of minerals within the top 50
inds_best <- order(res_CaMg)
inds_best50 <- which(inds_best < 51)
## Extract all phase names from the best 50
AllPhases <- unlist(lapply(filtered[inds_best50], function(x) rownames(x$pphases)))
## Count
table(AllPhases)

viz <- rbind(init, dat[which.min(res_CaMg),], cvec)

## making combinatorial plots for the best matched elemental solutions
par(mfrow=c(1,1))
out <- barplot(viz, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10),
               col=c("orange", "light green", "grey"))

legend("topleft", c("initial","rmse", "target"),
       fill=c("orange", "light green", "grey"), bty="n")

textCa <- filtered[[which.min(res_CaMg)]]$pphases
barplot(textCa$delta, names.arg=row.names(textCa), las=1, ylab="mol/kgw")
legend("bottomleft", c("Ca"))

################################################################

## Redox optimization

################################################################

ToOptimizePE <- function(pe, input, target) {
  tmp <- Distribute(input, prop="pe", values=pe)
  tmp <- AddProp(tmp, name="pe_Fix", values=paste(-pe, " O2(g)  0.5"), cat="pphases") 
  
  ## call PHREEQC and put the results into a list
  res <- .runPQC(tmp)
  
  ## we extract the total concentrations and pH and put them in a
  ## vector called "final"
  final <- ExtractComponents(res)
  
  ## make sure the "target" vector does not contain ()
  names(target) <- gsub("\\(.*$","", names(target))
  
  ## This function is written so that we can exclude some components
  ## from "target" by just removing them. However here we need to
  ## make sure that we compare the same components - compute the
  ## "intersection" of the vector names!
  components <- intersect(names(target), names(final))
  
  ## We compute the single numeric metric value. NB: I hard coded
  ## "rrmse" here but it can be changed
  ans <- rrmse(target[components], final[components], na.rm=FALSE)
  return(ans)
}

## We now use the R standard function "optimize()". Allowed range for
## pe is set to c(-8, 10)
#calib <- optimise(f=ToOptimizePE, input = NewInput, target = cvec,
#                  interval=c(-8,10), maximum = FALSE)## are the numbers subjective?
##Here we use the initial solution as the input
calib <- optimise(f=ToOptimizePE, input=inps[[1]], target = cvec,
                  interval=c(-7,7), maximum = FALSE)

## Look at the results
calib
str(calib)

temp <- Distribute(inps[[1]], prop="pe", values=calib$minimum) %>%
  AddProp(name="pe_Fix", values=paste(-calib$minimum, " O2(g)  0.5"), cat="pphases") 

## our final solution:
## require(magrittr) ## to use the %>% operator
## Calibrated <- Distribute(NewInput, prop="pe", values=calib$minimum) %>% .runPQC #this is originally commented out
#Calibrated <- .runPQC(Distribute(NewInput, prop="pe", values=calib$minimum))
Calibrated <- .runPQC(temp)
#InputFromList(Calibrated) this collects the parameter concentrations into a phreeqc script form

## comparison of the script before calibration and after

cviz <- rbind(ExtractComponents(Calibrated), ExtractComponents(filtered[[ind]]), target=cvec)
out <- barplot(cviz, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10),
               col=c("light green", "grey", "yellow"), las=1)
legend("topleft", c("Calibrated", "Final", "target"),
       fill=c("light green", "grey", "yellow"), bty="n")

######################################################################################

## another exemplary implementation varying Pyrite SI. This is applicable when there is a
## redox sensitive mineral in the mineral pool.

#####################################################################################

ToOptimizePyriteSI <- function(sipy, input, target) {
  tmp <- Distribute(input, prop="Pyrite", values=paste(sipy, "2"))
  
  ## call PHREEQC and put the results into a list
  res <- .runPQC(tmp)
  
  ## we extract the total concentrations and pH and put them in a
  ## vector called "final"
  final <- ExtractComponents(res)
  
  ## make sure the "target" vector does not contain ()
  names(target) <- gsub("\\(.*$","", names(target))
  
  ## This function is written so that we can exclude some components
  ## from "target" by just removing them. However here we need to
  ## make sure that we compare the same components - compute the
  ## "intersection" of the vector names!
  components <- intersect(names(target), names(final))
  
  ## We compute the single numeric metric value. NB: I hardcoded
  ## "rrmse" here but it can be changed
  ans <- rrmse(target[components], final[components], na.rm=FALSE)
  return(ans)
}

## We now use the R standard function "optimise()". Allowed range for
## pe is set to c(-8, 10)
calibpyr <- optimise(f=ToOptimizePyriteSI, input = NewInput, target = cvec,
                     interval=c(-8,8), maximum = FALSE)

## Look at the results
calibpyr
str(calibpyr)

## our final solution:
## require(magrittr) ## to use the %>% operator
## Calibrated <- Distribute(NewInput, prop="pe", values=calib$minimum) %>% .runPQC
CalibPyr <- .runPQC(Distribute(NewInput, prop="Pyrite", values=paste(calibpyr$minimum, "2")))

cviz2 <- rbind(ExtractComponents(CalibPyr), ExtractComponents(filtered[[ind]]), cvec)
out <- barplot(cviz2, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10),
               col=c("light green", "grey"), las=1)
