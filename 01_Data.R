
library(dplyr)
library(data.table)

## Baseline info ======================================================================================================
# 
DF_initial<-fread(file="/home/couchr02/Mendel_UKB/Source/Phenotype/March_2022_Update/ukb51070.tab", 
                  header = T, stringsAsFactors = F,
                  select = c("f.eid","f.53.0.0","f.21000.0.0","f.31.0.0","f.21022.0.0",
                             "f.30710.0.0","f.30870.0.0","f.30780.0.0","f.30760.0.0","f.30750.0.0",
                             "f.130706.0.0","f.130708.0.0","f.40001.0.0","f.40000.0.0","f.20116.0.0",
                             "f.189.0.0",
                             "f.93.0.0","f.94.0.0",
                             "f.4080.0.0","f.4079.0.0",
                             "f.21001.0.0","f.23104.0.0","f.50.0.0","f.21002.0.0",
                             "f.6177.0.0","f.6177.0.1","f.6177.0.2",
                             "f.6153.0.0","f.6153.0.1","f.6153.0.2",
                             "f.874.0.0","f.894.0.0","f.914.0.0",
                             "f.864.0.0","f.884.0.0","f.904.0.0",
                             "f.1289.0.0","f.1299.0.0","f.1309.0.0","f.1319.0.0",
                             "f.1329.0.0","f.1339.0.0","f.1349.0.0","f.1359.0.0",
                             "f.48.0.0","f.49.0.0"),
                  col.names = c("IID","REC_DATE","ETH","SEX","AGE",
                                "CRP","TG","LDL_C","HDL_C","HbA1C",
                                "D_T1D","D_T2D","CAUSE_Death","DATE_Death","Smoking",
                                "Townsend",
                                "SBP_M","DBP_M",
                                "SBP_A","DBP_A",
                                "BMI_M","BMI_B","height","weight",
                                "Rx_M_1","Rx_M_2","Rx_M_3",
                                "Rx_F_1","Rx_F_2","Rx_F_3",
                                "T_W","T_M","T_V",
                                "F_W","F_M","F_V",
                                "Cooked_vege","Raw_vege","Fresh_fruit","Dried_fruit",
                                "Oily_fish","Non_oily_fish","Processed_meat","Poultry",
                                "Waist_C","Hip_C"))

# Outcome ASCVD ------------------------------------------------------------------------------------------------------
DF_OUT<-fread(file="/home/couchr02/Mendel_Commun/Hasanga/UKB_Pheno/UKB_Data_MAJ_502K_07062022.txt", 
              header = T, stringsAsFactors = F,
              select = c("IID","DATE_ALL_MACE"))

DF_initial<-merge(x = DF_initial, y = DF_OUT, by = c("IID"))

# Sleep info ----------------------------------------------------------------------------------------------------------
DF_sleep<-fread(file="/home/couchr02/Mendel_UKB/Source/Phenotype/March_2022_Update/ukb51070.tab", 
                header = T, stringsAsFactors = F,
                select = c("f.eid","f.1160.0.0","f.1180.0.0","f.1200.0.0","f.1210.0.0"),
                col.names = c("IID","Sleep_duration","Chronotype","Insomnia","Snoring"))

# Dx sleep apnea
DF_Dx<-fread(file="/home/couchr02/Mendel_UKB/Source/Phenotype/March_2022_Update/ukb51070.tab", 
             header = T, stringsAsFactors = F,
             select = c("f.eid",paste0("f.41270.0.",0:242)),
             col.names = c("IID",paste0("f.41270.0.",0:242)))

DF_Dx$Sleep_apnea<-ifelse(grepl("G473",do.call(paste0,DF_Dx)),1,0)
# 10,108 participants 

DF_sleep<-merge(x = DF_sleep, y = DF_Dx[,c("IID","Sleep_apnea")], by = c("IID"))

DF_initial<-merge(x = DF_initial, y = DF_sleep, by = c("IID"))

# BP ------------------------------------------------------------------------------------------------------------------
DF_initial$SBP<-if_else(is.na(DF_initial$SBP_A),DF_initial$SBP_M,DF_initial$SBP_A)
DF_initial$DBP<-if_else(is.na(DF_initial$DBP_A),DF_initial$DBP_M,DF_initial$DBP_A)

DF_initial<-subset(DF_initial,select = -c(SBP_M,SBP_A,DBP_M,DBP_A))



# BMI -----------------------------------------------------------------------------------------------------------------
DF_initial$BMI<-if_else(is.na(DF_initial$BMI_M),DF_initial$BMI_B,DF_initial$BMI_M)

DF_initial<-subset(DF_initial,select = -c(BMI_M,BMI_B,height,weight))

# Rx ------------------------------------------------------------------------------------------------------------------

DF_initial$Rx_C<-if_else(DF_initial$Rx_M_1==1|DF_initial$Rx_M_2==1|DF_initial$Rx_M_3==1|DF_initial$Rx_F_1==1|DF_initial$Rx_F_2==1|DF_initial$Rx_F_3==1,1,0)
DF_initial$Rx_C<-if_else(is.na(DF_initial$Rx_C)&(DF_initial$Rx_F_1==-7|DF_initial$Rx_M_1==-7),0,DF_initial$Rx_C)
DF_initial$Rx_BP<-if_else(DF_initial$Rx_M_1==2|DF_initial$Rx_M_2==2|DF_initial$Rx_M_3==2|DF_initial$Rx_F_1==2|DF_initial$Rx_F_2==2|DF_initial$Rx_F_3==2,1,0)
DF_initial$Rx_BP<-if_else(is.na(DF_initial$Rx_BP)&(DF_initial$Rx_F_1==-7|DF_initial$Rx_M_1==-7),0,DF_initial$Rx_BP)
DF_initial$Rx_D<-if_else(DF_initial$Rx_M_1==3|DF_initial$Rx_M_2==3|DF_initial$Rx_M_3==3|DF_initial$Rx_F_1==3|DF_initial$Rx_F_2==3|DF_initial$Rx_F_3==3,1,0)
DF_initial$Rx_D<-if_else(is.na(DF_initial$Rx_D)&(DF_initial$Rx_F_1==-7|DF_initial$Rx_M_1==-7),0,DF_initial$Rx_D)
DF_initial$Rx_N<-if_else(is.na(DF_initial$Rx_M_1)&is.na(DF_initial$Rx_M_2)&is.na(DF_initial$Rx_M_3)&is.na(DF_initial$Rx_F_1)&is.na(DF_initial$Rx_F_2)&is.na(DF_initial$Rx_F_3),-3,1)
DF_initial$Rx_C<-if_else(is.na(DF_initial$Rx_C),0,DF_initial$Rx_C)
DF_initial$Rx_BP<-if_else(is.na(DF_initial$Rx_BP),0,DF_initial$Rx_BP)
DF_initial$Rx_D<-if_else(is.na(DF_initial$Rx_D),0,DF_initial$Rx_D)

DF_initial<-subset(DF_initial,select = -c(Rx_M_1,Rx_M_2,Rx_M_3,Rx_F_1,Rx_F_2,Rx_F_3))

save(DF_initial, file = "/mnt/sda/pauaud01/UKB_CVHS_ASCVD/Data/DF_initial.RData")

















