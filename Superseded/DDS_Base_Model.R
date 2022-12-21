setwd("C:/SANDEEP/SocioHydro_Model")

library(dataRetrieval)
library(Metrics)
library(corrplot)
library(tidyverse)
library(Hmisc)
source("SocHydModel_RainSurge.R")
# Returns numIter length list of entries to be peturbed
probPeturb<-function(x, numIter)
{
    # Input is xBounds & numIter.  
    # Returns numIter entry list with the indices which will be peturbed
    xDims<-nrow(x)
    probabilityVector<-1-log(1:numIter)/log(numIter)
    peturbIdx<-apply(matrix(unlist(lapply(probabilityVector, function(x) as.logical(rbinom(xDims, 1, x)))), byrow=TRUE, ncol=xDims), 1, which)
    return(peturbIdx)
}

#Read in index

fips = read.table("FIPS.txt")
Battery = read.csv("Max_Water_Level_Battery.csv")
precip = read.table("C:/SANDEEP/Precip/precip_timeseries.csv")

#DDS parameters
r= 0.3
numIter = 100000

#running code for census tracts of interest
for (j in 305:305)#length(fips$V1))
{
    #Read in index
    #Reading of fema claims is divided for HRE(295 census tracts) and LIS as the LIS filename and nrows is different than HRE
    if(j < 296)
    {fema_claims = read.table(paste("C:/SANDEEP/ForSandeep/FEMA_Claims_/HRE_and_LIS/FEMA_Claims_",fips$V1[j],".csv",sep = ''),sep = ',',skip = 1)
    fema_policies = read.table(paste("C:/SANDEEP/ForSandeep/FEMA_Policies_/HRE_and_LIS/HRE_FEMA_Policies_",fips$V1[j],".csv",sep = ''),sep = ',',skip = 1)}
    else
    {fema_claims = read.table(paste("C:/SANDEEP/ForSandeep/FEMA_Claims_/HRE_and_LIS/FEMAClaims_",fips$V1[j],".csv",sep = ''),sep = ',',skip = 2)
    fema_policies = read.table(paste("C:/SANDEEP/ForSandeep/FEMA_Policies_/HRE_and_LIS/FEMA_Policies_",fips$V1[j],".csv",sep = ''),sep = ',',skip = 1)}
    #Policy is missing one census tract(34017980100), a decoy value is added now should be removed later
    nhouses <- read.table(paste("C:/SANDEEP/Housing units&prices/Housing units/Sorted_by_year/final_processed/nhousing_",fips$V1[j],".csv",sep = ''),sep=',',skip = 1)
    vhouses <- read.table(paste("C:/SANDEEP/Housing units&prices/Housing prices/Sorted_by_year/final_processed/vhousing_",fips$V1[j],".csv",sep = ''),sep=',',skip = 1)
    #0 values of nhouses and vhouses are assigned NA
    nhouses$V3[nhouses$V3==0 | nhouses$V3=="null" | nhouses$V3=="-"] <- NA
    vhouses$V3[vhouses$V3==0 | vhouses$V3=="null" | vhouses$V3=="-"] <- NA
    if(is.na(mean(vhouses$V3)))
    {
        a=10000000 #there are about 20 tracts with no data so a random value is assigned for now
        Claims_Rescaled_v = fema_claims$V4/a
    }
    else 
    { Claims_Rescaled_v = fema_claims$V4/mean(vhouses$V3, na.rm = T)}
    
    if(is.na(mean(nhouses$V3)))
    {
        b=1000 #there are about 20 tracts with no data so a random value is assigned for now
        Claims_Rescaled_n = fema_claims$V3/b
        Claims_Rescaled_d = 0.1 * mean(nhouses$V3)/max(nhouses$V3) #Based on literature about 10 percent of total land in ct is used for residential
        Policy_Rescaled_a = fema_policies$V3/b
        
    }
    else 
    { Claims_Rescaled_n = fema_claims$V3/mean(nhouses$V3, na.rm = T)
    Claims_Rescaled_d = 0.1 * nhouses$V3/max(nhouses$V3)
    Policy_Rescaled_a = fema_policies$V3/mean(nhouses$V3, na.rm = T)
    }
    #Set calibration parameter boundaries
    xBounds.df = data.frame(matrix(ncol=2,nrow=11))
    colnames(xBounds.df)<-c("min", "max")
    
    #Surge threshold
    xBounds.df$min[1] = 0 #what about putting the min value of surge in lower bound(1)?
    xBounds.df$max[1] = 3.5 #wmax 
    
    #Rain threshold
    xBounds.df$min[2] = 0 #what about putting the min value of Rain in lower bound(0.03)?
    xBounds.df$max[2] = 0.15 #precipmax
    
    #alphad
    xBounds.df$min[3] = 0.001 
    xBounds.df$max[3] = 25 #constrained based on barplot from 100
    
    #alphaa
    xBounds.df$min[4] = 0.001
    xBounds.df$max[4] = 25
    
    #alphap
    xBounds.df$min[5] = 0.001
    xBounds.df$max[5] = 10 #increased to 10 from 3 based on sh output
    
    #alphar
    xBounds.df$min[6] = 0.001
    xBounds.df$max[6] = 1 #constrained based on barplot from 3
    
    #mewa
    xBounds.df$min[7] = 0.001
    xBounds.df$max[7] = 1 #constrained based on barplot from 3
    
    #mewp
    xBounds.df$min[8] = 0.001
    xBounds.df$max[8] = 2 #constrained based on barplot from 3
    
    #U_rate
    xBounds.df$min[9] = 0.01
    xBounds.df$max[9] = 3 #constrained based on barplot from 5
    
    #POT_S_max
    xBounds.df$min[10] = 0.001 
    xBounds.df$max[10] = 3.5 #max w - min threshold w
    
    #POT_R_max
    xBounds.df$min[11] = 0.001
    xBounds.df$max[11] = 0.15
    
    # Generate initial first guess
    x_init<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
    x_best = data.frame(x_init)
    RMSE_best <- 10000
        

    peturbIdx<-probPeturb(xBounds.df, numIter)
    # Peturb each entry by N(0,1)*r(x_max - x_min) reflecting if beyond boundaries
    sigma<-xBounds.df$max - xBounds.df$min
    
    for (i in 1:numIter)
    {
        # Set up test parameter values as x_test
        x_test<-as.matrix(x_best)
        
        # Get entries we will peturb
        idx<-peturbIdx[[i]]
        if (sum(idx) == 0) {idx = round(runif(1,1,11))}
        
        # Initialize vector of peturbations initially zeros with same length of x so we will add this vector to peturb x
        peturbVec<-rep(0, length(x_test))
        # Generate the required number of random normal variables
        N<-rnorm(length(x_test), mean=0, sd=1)
        
        # Set up vector of peturbations
        peturbVec[idx]<-r*N[idx]*sigma[idx]
        
        # Temporary resulting x value if we peturbed it
        testPeturb<-x_test + peturbVec  
        # Find the values in testPeturb that have boundary violations.  Store the indices in boundaryViolationsIdx
        boundaryViolationIdx<-which(testPeturb<xBounds.df$min | testPeturb > xBounds.df$max)
        
        # Reset those violated indices to the opposite peturbation direction
        peturbVec[boundaryViolationIdx]<-(-1*r*N[boundaryViolationIdx]*sigma[boundaryViolationIdx])
        
        # Find values still at violations of min or max and set them to the minimum or maximum values
        x_test<-x_test + peturbVec
        minViolationIdx<-which(x_test<xBounds.df$min)
        maxViolationIdx<-which(x_test>xBounds.df$max)
        x_test[minViolationIdx]<-xBounds.df$min[minViolationIdx]
        x_test[maxViolationIdx]<-xBounds.df$max[maxViolationIdx]
        
        #Socio-hydrological model    
        Results = SocHydModel_SurgePrecip(W = Battery$Max_Wl_Battery[10:52], Precip = precip$precip_max[10:52], x_test[1], x_test[2], x_test[3], x_test[4], x_test[5], x_test[6], 1, x_test[7], x_test[8], 
                                          1, 1, 1, 1, x_test[9], x_test[10], x_test[11])
        
        #Obj. Function 
        #RMSE_Test_N = rmse(Claims_Rescaled_n , Results$L ) / mean(Claims_Rescaled_n)
        RMSE_Test_V = rmse(Claims_Rescaled_v , Results$L ) / mean(Claims_Rescaled_v) #In our claim data, claim value is zero for several observations even when claim no is not, so it is not taken for calibration
        RMSE_Test_A = rmse(Policy_Rescaled_a , Results$A[31:43] ) / mean(Policy_Rescaled_a)
        RMSE_Test_D = rmse(Claims_Rescaled_d , Results$D[33:43] ) / mean(Claims_Rescaled_d)
        
        RMSE_Test <- sqrt(RMSE_Test_V^2 + RMSE_Test_A^2 + RMSE_Test_D^2 )
        #RMSE_Test <- RMSE_Test_V
        #Check if this simulation is better
        if (RMSE_Test < RMSE_best) 
        {
            x_best = x_test
            RMSE_best = RMSE_Test
            
        }
        
        
        #Print to console
        print_str = paste(fips$V1[j],"RMSE BEST ",i,":",RMSE_best)
        print(print_str)
        
    }
    
  
    #write results
    x_best <- round(x_best, digits = 5)
    colnames(x_best) <-"Best Value"
    rownames(x_best) <- c("Surge_threshold","Rain_threshold","alpha_d","alpha_a","alpha_p","alpha_r","mew_a","mew_p","U_rate","POT_S_max","POT_R_max")
    #write.table(x_best,paste("C:/SANDEEP/SocioHydro_Model/calibrated_parameters/SHparameters_",fips$V1[j],".csv",sep = ''),sep=",")
}

##Plot
par(mfrow=c(3,2))
#storm surge
plot(Battery$Year[10:52] ,Battery$Max_Wl_Battery[10:52], type = "o")


#awareness 
plot(Battery$Year[10:52] ,Results$A , type = "o", ylim = c(0,max(Policy_Rescaled_a)))
lines(Battery$Year[40:52],Policy_Rescaled_a, col = "red")


#precipitation
plot(precip$yr[10:52] ,precip$precip_max[10:52], type = "o" )

#preparedness
plot(Battery$Year[10:52] ,Results$P, type = "o" , ylim = c(0, max(Results$P)))

#loss
plot(Battery$Year[10:52] ,Results$L, type = "o")
lines(fema_claims$V2 , Claims_Rescaled_v, col = "red")

#density
plot(Battery$Year[10:52] ,Results$D, type = "o" , ylim = c(0,max(Claims_Rescaled_d)))
lines(Battery$Year[42:52], Claims_Rescaled_d, col = "red")
