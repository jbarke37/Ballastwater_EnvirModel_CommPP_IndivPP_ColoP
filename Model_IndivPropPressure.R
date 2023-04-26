### HOUSEKEEPING####

rm(list=ls())
set.seed(2935)
if (require(plyr) == FALSE){
install.packages("plyr")
}
library(plyr)

##### SIMULATION SET-UP PARAMETERS ########

ships_data <- read.csv("Discharge_volumes.csv")

##### SIMULATION SET-UP PARAMETERS ########
c=1
alpha_par1=0.005
alpha_par2=5

iterations=100
current_date=Sys.Date()

CPs <- c(1, 2, 5, 10, 30, 50, 100)  #colonization pressure
PPs <- c(1, 2, 5, 10, 30, 50, 100, 1000, 5000) #Community propagule pressure

# Create empty data frames to fill
risk_atleast <- data.frame()
risk_exactlyone <- data.frame()

##### SIMULATION ########

for (r in 1:nrow(ships_data))  { 
  for (nspec in CPs) {  #Community pressure, number of species
    for (PP in PPs)  { #number of organisms per species
      for (i in 1:iterations) {
        
        species_data=matrix(data=NA,nrow=as.numeric(nspec),ncol=11)
        colnames(species_data)=c("CP","PP","iter","c.val","spec.num","alpha.val","Comm.PP","indiv.PP","discharge_volume","P_Est","Est")
        species_data=as.data.frame(species_data)
        species_data$CP <- nspec
        species_data$PP <- PP
        species_data$iter <- i
        
        discharge_volume<-ships_data$DisVol[r]
        species_data$discharge_volume<-discharge_volume
        species_data$spec.num<-1:as.numeric(nspec)
        
        alpha_vals <- rbeta(nspec, shape1 = alpha_par1, shape2 = alpha_par2)
        
        species_data$alpha.val<-alpha_vals
        species_data$c.val<-c
        
        max_total_PP<- discharge_volume*PP #max number of organisms (volume * number of species per cubic meter of water )
        species_pool<-rep(unique(species_data$spec.num), each = floor((discharge_volume*PP)/nspec))  
        pool_summary=as.data.frame(table(species_pool))
        species_data$indiv.PP=pool_summary[,2]/discharge_volume
        species_data$Comm.PP <- sum(species_data$indiv.PP)
        species_data$P_Est=1-exp(1)^-(species_data$alpha.val*species_data$indiv.PP^species_data$c.val) #Calculate establishment probabilities
        species_data$Est=ifelse(species_data$P_Est<runif(nspec),0,1) # Probabilistically assign whether species establishes
        
        risk_atleast <- rbind(risk_atleast,
                              data.frame(Discharge_volume = discharge_volume,
                                         CommunityPP = PP,
                                         ColonizationP = nspec,
										 InivPP = PP/nspec,
                                         iter = i,
                                         Risk_Atleast = 1 - prod(1-species_data$P_Est)))
                              
          risk_exactlyone <- rbind(risk_exactlyone,
                              data.frame(Discharge_volume = discharge_volume,
                                    CommunityPP = PP,
                                     ColonizationP = nspec,
									 InivPP = PP/nspec,
                                    iter = i,
                                    n_SppEst = sum(species_data$Est)))
                                                       
           rm(species_data, pool_summary, potential_pool, max_total_PP, alpha_vals, discharge_volume)
                                                       
      }
    }
  }
}


write.csv(risk_atleast, "IndivPropPressure_Risk_AtleastOneSpp_Establishing.csv", row.names = FALSE)

if (require(dplyr) == FALSE){
  install.packages("dplyr")
}
library(dplyr)

risk_exactlyone_prob <- risk_exactlyone %>% 
  mutate(Group = case_when(n_SppEst == 0 ~ "None",
                           n_SppEst == 1 ~ "One",
                           n_SppEst > 1 ~ "More")) %>% 
  group_by(Discharge_volume, CommunityPP, ColonizationP, InivPP, Group) %>% 
  summarise(NSPPEST_Trip = n()) %>% 
  mutate(Prob = NSPPEST_Trip / iterations) %>% 
  filter(Group == "One")

write.csv(risk_exactlyone_prob , "IndivPropPressure_Risk_ExactlyOneSpp_Establishing.csv", row.names = FALSE)
