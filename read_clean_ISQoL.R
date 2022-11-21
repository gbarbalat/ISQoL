# read REHABASE
close.screen(all=TRUE)
rm(list=ls())
header=1;

library(dplyr)
library(tidyr)
library(mice)
library(parallel)
library(fastDummies)

##################
##################

#Load scores
First_Network_REHABASE_tot = read.csv(file="C:/Users/Guillaume/Downloads/extract_JULY_2022.csv",
                                      header=TRUE, sep=";",
                                      comment.char="#",encoding = "UTF-8")


mat_A_W_Y <- First_Network_REHABASE_tot %>%
  dplyr::select(Study.Subject.ID,CENTRE,
                Age, Sex,Education, #Employment,
                Dx, Dx2, Dx_SOMA, GAF, CGI, Fam, RQTH, #Parent, Adresseur, 
                Illness_Duration, First_Contact, N_Admissions, Marginalisation, Forensic, Addictions,
                SQoL18_SEL,SQoL18_RES,SQoL18_AUT,SQoL18_PHY,SQoL18_PSY,SQoL18_FRI,SQoL18_FAM,SQoL18_ROM,SQoL18_TOT,
                WEMWBS_TOT,
                IS_TOT,SERS_TOT,#IS_Sympt,IS_Disease,IS_Treatment,
                ISMI_TOT,#ISMI_Alien,ISMI_Approb,ISMI_Discrim,ISMI_Isol,ISMI_Intern,ISMI_Resist,
                STORI,
                #MAT, SIM, D2R_CCT, D2R_E, D2R_CC,
                # RLRI16_RLI, RLRI16_RLD, RLRI16_RTD, CVLT_RLI,CVLT_RLLT,MEMCHIF_MCD, MEMCHIF_MCI,
                # COMMISSIONS_time, COMMISSIONS_error,
                #ACSo_TOT,TIME
  ) %>%
  #mutate(across(c(Age,SQoL18_SEL:ACSo_TOT),~ gsub(",",".", .))) %>%
  mutate(across(c(Age,SQoL18_SEL:ISMI_TOT),~ as.numeric(.))) %>%
  rename_with(~gsub("SQOL18_","",.x)) %>%
  #filter(Study.Subject.ID=="01ABMA0280") %>%
  na_if("") %>%
  na_if(8888) #%>%
  ###
  # select(Study.Subject.ID,
  #        Diagnostic.psychiatrique.principal..codes.CIM10.reconnus.par.DSM5.,
  #        Dx) %>%
  ###
  # dplyr::group_by(Study.Subject.ID) %>%
  # unique()
  # dplyr::mutate(across(Age:ACSo_TOT,~first(na.omit(.)))) %>%
  # dplyr::mutate(across(Age:Addictions,~first(na.omit(.)))) %>%
  # dplyr::mutate(across(SQoL18_SEL:D2R_CC,~first(.))) %>%
  # dplyr::distinct(Study.Subject.ID,.keep_all = TRUE) 




##################
##################
#Clean A: binearize IS tot
mat_A_W_Y = mat_A_W_Y %>% 
  mutate(A=case_when(IS_TOT<=9 ~ as.character(0),
                     IS_TOT>9 ~ as.character(1))
         )

##################
##################
#Clean 
mat_A_W_Y<-mat_A_W_Y %>% 
  mutate(
    
    Sex=case_when(Sex=="Masculin"~"Male",
                  Sex=="F�minin"~"Female",
                  TRUE~NA_character_),
    Education=case_when(as.numeric(sapply(strsplit(Education,"=",fixed=TRUE),function(x) x[1]))>=12~">= 12 years",
                        Education=="Non Cotable/Non Pertinent"~NA_character_,
                        is.character(Education) ~ "< 12 years",
                        TRUE~NA_character_),
    # Employment=case_when(as.numeric(sapply(strsplit(Employment,",",fixed=TRUE),function(x) x[1]))==1~"EMPLOYED",
    #                      as.numeric(sapply(strsplit(Employment,",",fixed=TRUE),function(x) x[1]))==2~"EMPLOYED",
    #                      Employment=="99"~NA_character_,
    #                      is.character(Employment) ~ "UNEMPLOYED",
    #                      TRUE~NA_character_),
    N_Admissions=case_when(N_Admissions<=2~"<= 2",
                           N_Admissions>=3~">= 3",
                           #N_Admissions==2 | N_Admissions ==3~"2-3",
                           TRUE~NA_character_),
    First_Contact=case_when(First_Contact==">10 ans" ~ "> 10 years",
                            First_Contact=="5 � 10 ans" ~ "5 to 10 years",
                            First_Contact=="Non Cotable/Non Pertinent" ~NA_character_,
                            is.character(First_Contact)~"< 5 years",
                            TRUE~NA_character_),
    Fam=case_when(Fam=="Non demand� durant l'entretien" ~ NA_character_,
                  Fam=="Mari�(e)" ~ "In a relationship",
                  Fam=="PACSE" ~ "In a relationship",
                  Fam=="Union Libre" ~ "In a relationship",
                  is.character(Fam)~"Not in a relationship",
                  TRUE~NA_character_),
    Dx=case_when(Dx=="1-Troubles NEURODEVELOPPEMENTAUX" ~ "ASD",
                 Dx=="18-Troubles de la PERSONNALITE" ~ "PD",
                 Dx=="2-Spectre de la SCHIZOPHRENIE" ~ "SCZ",
                 Dx=="3-Troubles BIPOLAIRES" ~ "BAD",
                 Dx=="4-Troubles DEPRESSIFS" ~ "DEP",
                 Dx=="5-Troubles de l�ANXIETE" ~ "ANX" ,
                 Dx=="6-TOC" ~ "ANX",
                 Dx=="7-Troubles li�s � STRESS ou TRAUMATISME" ~ "ANX",
                 TRUE~as.character(Dx)),
    Dx_SOMA=case_when(Dx_SOMA=="Non demand� durant l'entretien" ~ NA_character_,
                      TRUE~as.character(Dx_SOMA)),
    Addictions=case_when(Addictions=="0" ~ "Addictions-",
                         is.character(Addictions) ~ "Addictions+",
                         TRUE~NA_character_),
    Forensic=case_when(Forensic=="Non" ~ "Forensic-",
                       Forensic=="Oui" ~ "Forensic+",
                       TRUE~NA_character_),
    Marginalisation=case_when(Marginalisation=="Non" ~ "Marginalisation-",
                              Marginalisation=="Oui" ~ "Marginalisation+",
                              TRUE~NA_character_),
    RQTH=case_when(RQTH=="Demande en attente" ~ "RQTH-",
                   RQTH=="Non" ~ "RQTH-",
                   RQTH=="Oui" ~ "RQTH+",
                   TRUE~NA_character_),
    # Adresseur=case_when(Adresseur=="Psychiatre libéral" ~ "Private",
    #                     Adresseur=="Psychologue libéral(e)" ~ "Private",
    #                     Adresseur=="Non demandé durant l'entretien" ~ NA_character_,
    #                     is.character(Adresseur) ~ "Public",
    #                     TRUE~NA_character_),
    STORI=case_when(STORI=="1 - MORATOIRE" ~ "MORATORIUM",
                    STORI=="2 - CONSCIENCE" ~ "REBUILDING" ,
                    STORI=="3 - PREPARATION" ~ "REBUILDING",
                    STORI=="4 - RECONSTRUCTION" ~ "REBUILDING",
                    STORI=="5 - CROISSANCE" ~ "GROWTH",
                    TRUE~NA_character_)
    
  )
#####################
#####################

mat_A_W=mat_A_W_Y

################################
################################
#make factor/ordinal factors
mat_A_W$Study.Subject.ID=factor(mat_A_W$Study.Subject.ID)
mat_A_W$CENTRE=factor(mat_A_W$CENTRE)

mat_A_W$Sex=factor(mat_A_W$Sex)
mat_A_W$Education=factor(mat_A_W$Education)
#mat_A_W$Employment=factor(mat_A_W$Employment)
mat_A_W$RQTH=factor(mat_A_W$RQTH)
mat_A_W$Fam=factor(mat_A_W$Fam)
#mat_A_W$Adresseur=factor(mat_A_W$Adresseur)
mat_A_W$Marginalisation=factor(mat_A_W$Marginalisation)
mat_A_W$Forensic=factor(mat_A_W$Forensic)
mat_A_W$Dx=factor(mat_A_W$Dx)
mat_A_W$Dx2=factor(mat_A_W$Dx2)
mat_A_W$Dx_SOMA=factor(mat_A_W$Dx_SOMA)
mat_A_W$Addictions=factor(mat_A_W$Addictions)


mat_A_W$N_Admissions <- factor(mat_A_W$N_Admissions, 
                                          order = FALSE, #TRUE  makes it an ordinal variable
                                          #levels = c("< 2","2-3","> 3"))
                                          levels = c("<= 2",">= 3"))

mat_A_W$First_Contact <- factor(mat_A_W$First_Contact, 
                                           order = FALSE, #this makes it an ordinal variable
                                           levels = c("< 5 years","5 to 10 years","> 10 years"))

mat_A_W$STORI <- factor(mat_A_W$STORI, 
                                   order = FALSE, #this makes it an ordinal variable
                                   levels = c("MORATORIUM",
                                              "REBUILDING",
                                              "GROWTH"))
#IS_Sympt
#IS_Disease

#no need to transform the below into factors
#WEMWBS
#SQoL
#ISMI
#IS_TOT
################################
################################



#######################
#######################
#make all variables as dummies
# mat_A_W_alldum_pre <- mat_A_W %>% 
#   ungroup() %>%
#   
#   dplyr::select(c(Study.Subject.ID,
#                   Age,CENTRE,
#                   Sex,Education,Employment,
#                   RQTH,Fam,Adresseur,Marginalisation,Forensic,
#                   Dx,Dx2,Dx_SOMA,Addictions,
#                   GAF,CGI,#N_PEC_bin,
#                   First_Contact,N_Admissions,
#                   
#                   SQoL18_SEL:SQoL18_TOT,
#                   WEMWBS_TOT,
#                   ISMI_TOT,
#                   IS_TOT,
#                   STORI)) %>%
#   mutate(Age = case_when(Age <=35 ~ as.character(0),
#                       Age >35 ~ as.character(1)),
#        GAF = case_when(GAF <=60 ~ as.character(0),
#                        GAF>60 ~ as.character(1)),
#        CGI = case_when(CGI<4 ~ as.character(0),
#                        CGI>=4 ~ as.character(1)),
#        across(c(SQoL18_SEL:SQoL18_TOT), ~ case_when(.x<=50 ~ as.character(0),
#                                                     .x>50 ~ as.character(1))),
#        WEMWBS_TOT =case_when(WEMWBS_TOT<=42 ~ as.character(0),
#                              WEMWBS_TOT>42 ~ as.character(1)),
#        ISMI_TOT = case_when(ISMI_TOT<=2.5 ~ as.character(0),
#                             ISMI_TOT>2.5 ~ as.character(1)),
#        IS_TOT= case_when(IS_TOT<=9 ~ as.character(0),
#                          IS_TOT>9 ~ as.character(1))
#        )


#make NA as an additional column (keeps more data)
# mat_A_W_alldum_NA <- dummy_cols(
#   .data=  dplyr::select(is.factor(mat_A_W)),
#   select_columns = NULL,
#   remove_first_dummy = TRUE,#to avoid multicol
#   remove_most_frequent_dummy = FALSE,
#   ignore_na = FALSE,
#   split = NULL,
#   remove_selected_columns = TRUE#removes the columns used to generate the dummy columns
#   ) %>%
#   mutate(across(everything(.),
#                 ~ replace(., is.na(.), median(.))
#                 )
#          )


################################
################################


##########
#SAVE
##########
save(mat_A_W,
     file="mat_A_W_ISQoL.RData")
##########