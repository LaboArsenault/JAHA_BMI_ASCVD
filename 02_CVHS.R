library("dplyr")
library("data.table")
library("stringr")


load("/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Data/DF_initial.RData")

# Classification eth (1878 Don't Know ou prefer not to answer) --------------------------------------------------------
DF_initial$ETH<-ifelse(DF_initial$ETH<0,0,DF_initial$ETH)
DF_initial$ETH<-str_extract(DF_initial$ETH,"\\d")


DF_All<-DF_initial

## Physical activity score --------------------------------------------------------------------------------------------
DF_All$T_W<-ifelse(is.na(DF_All$T_W),0,DF_All$T_W)
DF_All$T_M<-ifelse(is.na(DF_All$T_M),0,DF_All$T_M)
DF_All$T_V<-ifelse(is.na(DF_All$T_V),0,DF_All$T_V)
DF_All$F_W<-ifelse(is.na(DF_All$F_W),0,DF_All$F_W)
DF_All$F_M<-ifelse(is.na(DF_All$F_M),0,DF_All$F_M)
DF_All$F_V<-ifelse(is.na(DF_All$F_V),0,DF_All$F_V)

DF_All$T_W<-ifelse((DF_All$T_W==-1|DF_All$T_W==-3)&(DF_All$F_W>0),10,DF_All$T_W)
DF_All$T_M<-ifelse((DF_All$T_M==-1|DF_All$T_M==-3)&(DF_All$F_M>0),10,DF_All$T_M)
DF_All$T_V<-ifelse((DF_All$T_V==-1|DF_All$T_V==-3)&(DF_All$F_V>0),10,DF_All$T_V)

DF_All$F_W<-ifelse((DF_All$F_W==-1|DF_All$F_W==-3)&(DF_All$T_W==0),0,DF_All$F_W)
DF_All$F_M<-ifelse((DF_All$F_M==-1|DF_All$F_M==-3)&(DF_All$T_M==0),0,DF_All$F_M)
DF_All$F_V<-ifelse((DF_All$F_V==-1|DF_All$F_V==-3)&(DF_All$T_V==0),0,DF_All$F_V)


# Temps total d'AP
DF_All$Tot<-rowSums(data.frame(DF_All$T_W,DF_All$T_M,DF_All$T_V))
DF_All$Tot_C<-if_else(DF_All$Tot>=960,0,1)
DF_All<-DF_All%>%dplyr::filter(Tot_C==1)
# - 1814 participants

DF_All$T_W<-ifelse(DF_All$T_W<10,0,DF_All$T_W)
DF_All$T_M<-ifelse(DF_All$T_M<10,0,DF_All$T_M)
DF_All$T_V<-ifelse(DF_All$T_V<10,0,DF_All$T_V)
DF_All$T_W<-ifelse(DF_All$T_W>=180,180,DF_All$T_W)
DF_All$T_M<-ifelse(DF_All$T_M>=180,180,DF_All$T_M)
DF_All$T_V<-ifelse(DF_All$T_V>=180,180,DF_All$T_V)



# Calcul MET min 
DF_All$MET_W<-(3.3*DF_All$T_W*DF_All$F_W)
DF_All$MET_M<-(4*DF_All$T_M*DF_All$F_M)
DF_All$MET_V<-(8*DF_All$T_V*DF_All$F_V)
DF_All$MET<-rowSums(data.frame(DF_All$MET_W,DF_All$MET_M,DF_All$MET_V)) #ALL
DF_All$MET_M_V<-rowSums(data.frame(DF_All$MET_M,DF_All$MET_V)) #Mod et Vig seulement

# Score physical activity
DF_All$S_AP<-if_else(DF_All$MET_M_V>=600,0,1)

DF_All<-subset(DF_All, select = -c(T_W,T_M,T_V,F_W,F_M,F_V,Tot,Tot_C,MET_W,MET_M,MET_V,MET,MET_M_V))

# Waist-to-hip ratio --------------------------------------------------------------------------------------------------
DF_All$WHR<-DF_All$Waist_C/DF_All$Hip_C
DF_All<-subset(DF_All, select = -c(Waist_C,Hip_C))



# Remove NA -----------------------------------------------------------------------------------------------------------
DF_All<-DF_All%>%dplyr::filter(!is.na(SBP)&!is.na(DBP)&!is.na(CRP)&
                                 !is.na(TG)&!is.na(LDL_C)&!is.na(HDL_C)&
                                 !is.na(HbA1C)&!is.na(BMI)&ETH>0&!is.na(WHR))

DF_All$T1D<-if_else(!is.na(DF_All$D_T1D)&is.na(DF_All$D_T2D),1,0)

##### Metabolic score #####
# Score BP ------------------------------------------------------------------------------------------------------------
DF_All$S_BP<-if_else(DF_All$SBP<130&DF_All$DBP<=80&DF_All$Rx_BP==0,0,1)

# Score CRP -----------------------------------------------------------------------------------------------------------
DF_All$S_CRP<-if_else(DF_All$CRP<3|DF_All$CRP>200,0,1) 

# Score TG ------------------------------------------------------------------------------------------------------------
DF_All$S_TG<-if_else(DF_All$TG<2.3,0,1)

# Score LDL -----------------------------------------------------------------------------------------------------------
DF_All$S_LDL_C<-if_else(DF_All$LDL_C<3&DF_All$Rx_C==0,0,1)

# Score HDL -----------------------------------------------------------------------------------------------------------
DF_All$S_HDL_C<-if_else(DF_All$HDL_C>1,0,1)

# Score HbA1C ---------------------------------------------------------------------------------------------------------
DF_All$S_HbA1C<-if_else(DF_All$HbA1C<42&DF_All$Rx_D==0,0,1)
DF_All$S_HbA1C<-if_else(DF_All$S_HbA1C==1&DF_All$T1D==1,0,DF_All$S_HbA1C)

# Score total ---------------------------------------------------------------------------------------------------------
DF_All$Score_metabolique<-rowSums(data.frame(DF_All$S_BP,DF_All$S_CRP,DF_All$S_TG,DF_All$S_LDL_C,DF_All$S_HDL_C,
                                             DF_All$S_HbA1C))
DF_All$Score_metabolique_factor<-as.factor(DF_All$Score_metabolique)


# Score Smoking -------------------------------------------------------------------------------------------------------
DF_All<-dplyr::filter(DF_All,Smoking>=0) #Aucun NA donc filtre pour -3 (ne veux pas répondre)
DF_All$Smoking<-if_else(DF_All$Smoking==0|DF_All$Smoking==1,0,1)

DF_All<-DF_All%>%dplyr::filter(Processed_meat>=0,Poultry>=0,Oily_fish>=0,Non_oily_fish>=0,
                               Cooked_vege>=0,Raw_vege>=0,Fresh_fruit>=0,Dried_fruit>=0)

DF_All$Sums_F_L<-rowSums(data.frame(DF_All$Cooked_vege, DF_All$Raw_vege, DF_All$Dried_fruit/2, DF_All$Fresh_fruit))
DF_All$S_F_L<-if_else(DF_All$Sums_F_L<5,1,0)

# Sleep score ---------------------------------------------------------------------------------------------------------
DF_All<-DF_All%>%filter(Sleep_duration>0&Insomnia>0)

DF_All$C_Sleep_duration<-if_else(DF_All$Sleep_duration==7|DF_All$Sleep_duration==8|DF_All$Sleep_duration==9,0,1)
DF_All$C_Insomnia<-if_else(DF_All$Insomnia==1|DF_All$Insomnia==2,0,1)

DF_All$Sleep_Score<-DF_All$C_Sleep_duration+DF_All$C_Insomnia+DF_All$Sleep_apnea
DF_All$Sleep_Score_2<-if_else(DF_All$Sleep_Score==0,0,1)

# CVHS ----------------------------------------------------------------------------------------------------------------
DF_All$CVHS<-rowSums(data.frame(DF_All$Score_metabolique, DF_All$S_AP, DF_All$S_F_L, DF_All$Sleep_Score_2, DF_All$Smoking))


DF_All<-subset(DF_All, select = -c(Cooked_vege, Raw_vege, Dried_fruit, Fresh_fruit, Sleep_duration, Insomnia, C_Sleep_duration,
                                   C_Insomnia, Sleep_apnea, Chronotype, Snoring, D_T1D, Sums_F_L))

save(DF_All, file = "/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Data/DF_All_02.RData")
