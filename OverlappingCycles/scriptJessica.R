library(BreedingSchemeLanguage)
library(reshape)

##Parameter setting
#parameters that are under breeder control
nparents<- 30 
ncand<- 1000 #number tested in the 1st year
ncyc<- 15 #keep same, number of cycles
nself<- 4 # number of selfing generations cycles
advance1<- 200 #number tested in the 2nd year
advance2<- 100 #number tested in the 3rd year
loc1<-1 
loc2<-c(1:4)
loc3<-c(1:4)
rep1<- 1
rep2<- 2
rep3<- 2
plottype1<- 'Reduced'
plottype2<- 'Standard'
plottype3<- 'Standard'
#nchecks<- 5 #number tested each year as checks

#parameters for advanced users
NeSpecies<- 980
NePed<- 33
EqCg<- 5
nQTL<- 500
npbFounders<- 140

#calculated parameters
cycDur<- nself+3 #number of selfing+testing seasons
nbase<- cycDur*nparents

#Simluate the founders of IRRI germplasm
simEnv <- defineSpecies(nSim = 1, nChr = 12, effPopSize = NeSpecies,
                       nQTL=nQTL,lengthChr=110,
                       saveDataFileName = 'Indica500qtl')
# simEnv <- defineSpecies(loadData="Indica500qtl")
defineVariances(gVariance=1,gByLocVar = 0.75, 
                gByYearVar = 0.75, 
                fracGxEAdd = 1,
                plotTypeErrVars = c(Standard = 0.5, Reduced=1))
initializePopulation(sEnv = simEnv, nInd = npbFounders) #pop 0
popno<-0

##Simulate current IRRI breeding materials
#Wright-fisher population to create population-wide inbreeding
#assumes selection was not effective, so drift is the main force
for(i in 1:EqCg){ 
  cross(nProgeny=NePed) #Effective population size of 33 for ~5 gen.
  popno<- popno+1
}
cross(nProgeny=nbase) 
popno<- popno+1

#Create fixed lines to be sampled from to start the program
for(i in 1:nself){
  selfFertilize(sEnv = simEnv, nProgeny = nbase, popID = NULL) #S1
  popno<- popno+1
}
popBase<- popno

##Start of cyclical breeding process
popnos<- c()
popcand<-c()
activ<-c()
seasons<-c()
for(i in 1:ncyc){#start the loop
  if(i<cycDur){ #sample parents from base until tested lines are avaliable
    BreedingSchemeLanguage::select(nSelect=nparents, # random pick from base  
                                   popID=popBase, random=T) 
  }else{
    BreedingSchemeLanguage::select(nSelect=nparents, popID=popcand) 
  }
  popno<- popno+1  
  cross(nProgeny=ncand, popID=popno) #starting population
  popno<- popno+1  
  popnos<- c(popnos, popno)
  popnoX<- popno #popno of the cross
  season<- i #crossing season
  seasons<- append(seasons, season)
  activ<- append(activ,'crossing')

  #selfing generations
  for(k in 1:nself){
    if(k==1){
      selfFertilize(sEnv = simEnv, nProgeny = ncand, popID = popnoX)
    }else{
      selfFertilize(sEnv = simEnv, nProgeny = ncand, popID = popno) 
    }
  popno<- popno+1
  activ<- append(activ, 'selfing')
  season<- season+1
  seasons<- append(seasons, season)
  }

  #First testing stage
  season<- season+1
  seasons<- append(seasons, season)
  activ<- append(activ, 'test1')
  phenotype(nRep=rep1, plotType=plottype1, locations=loc1, years=season) #season 6 is the first 
  BreedingSchemeLanguage::select(nSelect=advance1, popID=popno) # select based on pheno
  popno<- popno+1 
  popnos<- c(popnos, popno)
  popcand<- append(popcand, popno) 

  #Second stesting stage
  season<- season+1
  seasons<- append(seasons, season)
  activ<- append(activ, 'test2')
  phenotype(nRep=rep2, plotType=plottype2, locations=loc2, years=season) #advanced lines 1
  predictValue(popID=popno) #blup for selection among current candidates
  BreedingSchemeLanguage::select(nSelect=advance2, popID=popno) 
  popno<- popno+1
  popnos<- c(popnos, popno)
  popcand<- append(popcand, popno) 

  #Third testing stage
  season<- season+1
  seasons<- append(seasons, season)
  activ<- append(activ, 'test3')
  phenotype(nRep=rep3, plotType=plottype3, locations=loc3, years=season) #advanced lines 2
  predictValue(popID=popcand) #blup for selection among ALL candidates
  BreedingSchemeLanguage::select(nSelect=nparents, popID=popcand)
  popno<- popno+1
  popnos<- c(popnos, popno)
  popparents<- popno
}#End of breeding program

#Summary of activities and pops by season
tab<- data.frame(seasons, activ, popnos)
cast(activ~seasons, value='popnos', data=tab)

outputResults(summarize=F,saveDataFileName='test')



