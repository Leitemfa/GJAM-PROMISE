#Interpreting the factorial analysis
#Ref: http://genomicsclass.github.io/book/pages/interactions_and_contrasts.html

GJAM.Factorial <- function(output,F1Name,F2Name,Taxa,ng,burning) {
  # GJAM.Factorial(output_all.F,F1Name = "Inoculation",F2Name = "Variety", 
  #                Taxa=ExName, ng=2000,burning=500)
  #Interpreting the two-way factorial design from GJAM
  # Fact1Name <- fact1
  # Fact2Name <- fact2
  require(gjam)
  require(tidyr)
  require(dplyr)
  Fact1Name <- F1Name
  Fact2Name <- F2Name
  outputGJAM <- output
  # outputGJAM <- selected.model
  TaxaName <- Taxa
  # TaxaName <- ExName
  # ng <- 2000
  # burning <- ng-100
  
  # print(colnames(outputGJAM$prediction$ypredMu)
  xdata.Fac <- as.data.frame(outputGJAM$inputs$xdata)
  colnames(xdata.Fac) <- colnames(outputGJAM$inputs$xdata)
  Fact.Treat <- xdata.Fac[,c(Fact1Name ,Fact2Name )]
  fact1 <- levels(xdata.Fac[,Fact1Name])
  fact2 <- levels(xdata.Fac[,Fact2Name])
  
  
  # Fact.Treat.Vector <- paste(levels(output$inputs$xdata[,3]),levels(output$inputs$xdata[,5]))
  Fact.Treat.Code <- data.frame(Fact1 = rep(fact1,length(fact2)),Fact2 = rep(fact2,each = length(fact1))) %>% 
    unite("TreatCode",Fact1:Fact2)
  
  TaxaSelected <- paste0(TaxaName,"_")
  df.gibbsT <- outputGJAM$chains$bgibbs[burning:ng,] %>% 
    data.frame()
  colnames(df.gibbsT) <- colnames(outputGJAM$chains$bgibbs[burning:ng,])
  df.gibbs <- df.gibbsT %>% 
    dplyr::select(starts_with(TaxaSelected)) #%>% colnames
  Interaction.Terms <- grep("\\:",colnames(df.gibbs))
  MainEffects <- colnames(df.gibbs)[-Interaction.Terms]
  EffectsF1 <- paste0(TaxaSelected,Fact1Name,unique(fact1))
  EffectsF2 <- paste0(TaxaSelected,Fact2Name,unique(fact2))
  
  # colnames(df.gibbs)[MainEffects]
  # df.gibbs[,EffectsF1[-1]] <- df.gibbs[,EffectsF1[-1]] + df.gibbs[,1]
  # boxplot(df.gibbs[,MainEffects])
  # boxplot(output$chains$bgibbs[,1:8])
  modelformula <- model.matrix(as.formula(paste("~", Fact1Name, "*",Fact2Name)),xdata.Fac)

  #simple way
  dd <- xdata.Fac
  UniqueMat <- unique(modelformula)
  rownames(UniqueMat) <- unique(paste(paste0(Fact1Name,dd[,Fact1Name]),paste0(Fact2Name,dd[,Fact2Name]),sep = ":"))
  AdjMat <- t(UniqueMat)
  
    # AdjMat <- t(modelformula) %*% modelformula
  # AdjMat[AdjMat>0] <- 1
  # AdjMat[lower.tri(AdjMat)] <- 0
  
  # colnames(df.gibbs)[Interaction.Terms]
  # df.gibbs[,Interaction.Terms] <- df.gibbs[,EffectsF1[-1]]+
  #   (df.gibbs[,colnames(df.gibbs)[Interaction.Terms]]+df.gibbs[,EffectsF2[-1]])
  
  #Build a matrix of recontrast
  df.gibbs2 <- data.frame(as.matrix(df.gibbs) %*% AdjMat)
  # boxplot(df.gibbs[,colnames(df.gibbs)[Interaction.Terms]])
  # boxplot(df.gibbs[,colnames(df.gibbs)])
  # 
  # boxplot(df.gibbs[,colnames(df.gibbs)[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16)]],las=2, cex.axis=0.8)
  # abline(h = 0)
  
  
  colnames(df.gibbs2) <- unique(paste(dd[,Fact1Name],dd[,Fact2Name],sep = "_")) 
  return(df.gibbs2)
}

#Outputs Dorota
outputBactBK2week
outputBactEndo2week
outputBactRh2week

outputFunBK2week
outputFunEndo2week
outputFunRh2week

outputBactBK3week
outputBactEndo3week
outputBactEndoSand3week
outputBactRh3week

outputFunBK3week
outputFunEndo3week
outputFunRh3week


#Example
library(gjam)

selected.model <- outputBactEndoSand3week

dim(selected.model$chains$bgibbs)
ExName <- colnames(selected.model$parameters$betaMu)[1]
# paste0(TaxaSelected,Fact1Name,unique(fact1))
library(dplyr)
ng <- dim(selected.model$chains$bgibbs)[1]
fact1 <- "Soil"
fact2 <- "Striga"
df.gibbs3 <- GJAM.Factorial(selected.model,F1Name = fact1,F2Name = fact2, 
                            Taxa=ExName, ng=ng,burning=burnin)
colnames(df.gibbs3)

library(ggplot2)
df.gibbs3 %>% 
  gather(key = "VarName", value = "Coeff") %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  ggplot(aes(x = Fact1, y = Coeff, fill = Fact2)) + 
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  theme_bw()

df.gibbs3 %>% 
  gather(key = "VarName", value = "Coeff") %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  ggplot(aes(x = Fact1, y = Coeff, fill = Fact2)) + 
  geom_boxplot(outlier.size = 0) + 
  stat_summary(fun.y="mean", geom="point", shape=23, size=4, fill="white") +
  guides(fill=FALSE)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  theme_bw()

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

lvF1 <- levels(as.data.frame(selected.model$inputs$xdata)[,fact1])
lvF2 <- levels(as.data.frame(selected.model$inputs$xdata)[,fact2])

df.gibbs3 %>% 
  gather(key = "VarName", value = "Coeff") %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% 
  ggplot(aes(x = Fact2, y = Coeff, fill = Fact1)) +
  # guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
  stat_summary(fun.y="median", geom="point", shape=23, size=4, fill="white") +
  stat_smooth(method="lm", formula=y~1, se=FALSE, aes(group = 1)) +
  # stat_summary(fun.y = "mean", color = "red", geom = "line", aes(group = 1))+
  # geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  theme_bw()

df.gibbs3 %>% 
  gather(key = "VarName", value = "Coeff") %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% 
  ggplot(aes(x = Fact1, y = Coeff, fill = Fact2)) +
  # guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
  stat_summary(fun.y="median", geom="point", shape=23, size=4, fill="white") +
  stat_smooth(method="lm", formula=y~1, se=FALSE, aes(group = 1)) +
  # stat_summary(fun.y = "mean", color = "red", geom = "line", aes(group = 1))+
  # geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  theme_bw()

barplot(apply(df.gibbs3,2,mean))
barplot(output1$parameters$betaMu[,1])

multinames <- function(varnames){
  print(varnames)
  res <- GJAM.Factorial(output = selected.model,
                        F1Name = fact1,F2Name = fact2,
                        Taxa=varnames, ng=ng,burning=burnin)
  return(res)
}

Fact.Res <- lapply(colnames(selected.model$parameters$betaMu),multinames)
names(Fact.Res) <- colnames(selected.model$parameters$betaMu)
titleName <- names(Fact.Res[110])
Fact.Res[[110]]  %>% 
  gather(key = "VarName", value = "Coeff") %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% 
  ggplot(aes(x = Fact2, y = Coeff, fill = Fact1)) +
  # guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
  stat_summary(fun.y="median", geom="point", shape=23, size=4, fill="white") +
  stat_smooth(method="lm", formula=y~1, se=FALSE, aes(group = 1)) +
  # stat_summary(fun.y = "mean", color = "red", geom = "line", aes(group = 1))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # ggtitle(titleName)+
  # scale_x_discrete(limits =fact1)+
  theme_bw()

Fact.Res[[110]]  %>% 
  gather(key = "VarName", value = "Coeff") %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% 
  ggplot(aes(x = Fact2, y = Coeff, fill = Fact1)) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
  stat_summary(fun.y="median", geom="point", shape=23, size=4, fill="white") +
  stat_smooth(method="lm", formula=y~1, se=FALSE, aes(group = 1)) +
  # stat_summary(fun.y = "mean", color = "red", geom = "line", aes(group = 1))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  theme_bw()

#All Variables
do.call(rbind,Fact.Res[1:2])
Fact.Res.DF <- data.frame(TaxVar = rep(names(Fact.Res),each=nrow(data.frame(Fact.Res[1]))),
                          do.call(rbind,Fact.Res))
View(Fact.Res.DF)

#Selecionar taxas específicos
TaxVarS <- unique(Fact.Res.DF$TaxVar)
Fact.Res.DF  %>% 
  gather(key = "VarName", value = "Coeff", -TaxVar) %>%
  filter(TaxVar %in% TaxVarS[1:9]) %>% #analisa de 1 a 9
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% 
  ggplot(aes(x = Fact1, y = Coeff, fill = Fact2)) +
  guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
  stat_summary(fun.y="median", geom="point", shape=23, size=4, fill="white") +
  stat_smooth(method="lm", formula=y~1, se=FALSE, aes(group = 1)) +
  # stat_summary(fun.y = "mean", color = "red", geom = "line", aes(group = 1))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  facet_wrap(~TaxVar)+
  theme_bw()

# Fact.Res.DF.adj <- Fact.Res.DF
# Fact.Res.DF.adj[,-1] <- Fact.Res.DF.adj[,-1]-rowMeans(Fact.Res.DF.adj[,-1])
# t1 <- Fact.Res.DF.adj  %>% # pegar os significativos
#   gather(key = "VarName", value = "Coeff", -TaxVar) %>%View
#   group_by(TaxVar,VarName) %>% 
#   summarise(Median = median(Coeff),Q1 = quantile(Coeff, 0.025), 
#             Q3 = quantile(Coeff, 0.975)) %>% 
#   ungroup() %>% 
#   mutate(Sign = ifelse(Q1<0 & Q3>0, "", "*")) %>% 
#   filter(Sign == "*")
# 
# SigVars <- as.vector(unique(t1$TaxVar))
# boxplot(apply(Fact.Res.DF[,-1],1,mean))

Fact.Res.DF[,paste0(lvF1[1],"_",lvF2[1])]


tH <- Fact.Res.DF  %>% # estimate the coefficients
  gather(key = "VarName", value = "Coeff", -TaxVar) %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge", remove = FALSE) %>% #View
  group_by(TaxVar,Fact1,Fact2,VarName) %>% 
  summarise(Mean = mean(Coeff),Q1 = quantile(Coeff, 0.025), 
            Q3 = quantile(Coeff, 0.975)) %>% 
  ungroup() 
  
View(tH)

tH0 <- tH %>% 
  group_by(TaxVar) %>% 
  summarise(H0 = mean(Mean)) %>% 
  ungroup()

tHF1 <- tH %>% 
  group_by(TaxVar,Fact1) %>% 
  summarise(H0F1 = mean(Mean)) %>% 
  ungroup()

tHF2 <- tH %>% 
  group_by(TaxVar,Fact2) %>% 
  summarise(H0F2 = mean(Mean)) %>% 
  ungroup()

tHcF1 <- Fact.Res.DF  %>% # estimate the coefficients
  gather(key = "VarName", value = "Coeff", -TaxVar) %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge", remove = FALSE) %>% #View
  group_by(TaxVar,Fact1,Fact2,VarName) %>% 
  summarise(Mean = mean(Coeff),Q1 = quantile(Coeff, 0.025), 
            Q3 = quantile(Coeff, 0.975)) %>% 
  ungroup() %>% 
  filter(Fact1 == lvF1[1]) %>% 
  select(TaxVar,Fact2,Mean) %>% 
  rename(HcF1 = Mean)
tHcF2 <- Fact.Res.DF  %>% # estimate the coefficients
  gather(key = "VarName", value = "Coeff", -TaxVar) %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge", remove = FALSE) %>% #View
  group_by(TaxVar,Fact1,Fact2,VarName) %>% 
  summarise(Mean = mean(Coeff),Q1 = quantile(Coeff, 0.025), 
            Q3 = quantile(Coeff, 0.975)) %>% 
  ungroup() %>% 
  filter(Fact2 == lvF2[1]) %>% 
  select(TaxVar,Fact1,Mean) %>% 
  rename(HcF2 = Mean)

tHF <- tH %>% 
  left_join(tH0, by = "TaxVar") %>% 
  left_join(tHF1, by = c("TaxVar","Fact1")) %>% 
  left_join(tHF2, by = c("TaxVar","Fact2")) %>% 
  left_join(tHcF1, by = c("TaxVar","Fact2")) %>% 
  left_join(tHcF2, by = c("TaxVar","Fact1")) %>% 
  mutate(SignH0 = ifelse(Q1<=H0 & Q3>=H0, "", "*")) %>%  #If any difference between treatments
  mutate(SignH0F1 = ifelse(Q1<=H0F1 & Q3>=H0F1, "", "*")) %>% #if any difference between Fact1 treatments
  mutate(SignH0F2 = ifelse(Q1<=H0F2 & Q3>=H0F2, "", "*")) %>%  #if any difference between Fact2 treatments
  mutate(SignHcF1 = ifelse(Q1<=HcF1 & Q3>=HcF1, "", "*")) %>%  #if treatments in Fact1 differs from the control within Fact2
  mutate(SignHcF2 = ifelse(Q1<=HcF2 & Q3>=HcF2, "", "*")) #if treatments in Fact2 differs from the control within Fact1

View(tHF)
dim(tHF)

(SigVarsH0 <- as.vector(unique(as.vector(tHF[tHF$SignH0=="*","TaxVar"]))))
(SigVarsH0F1 <- as.vector(unique(as.vector(tHF[tHF$SignH0F1=="*","TaxVar"]))))
(SigVarsH0F2 <- as.vector(unique(as.vector(tHF[tHF$SignH0F2=="*","TaxVar"]))))
(SigVarsHcF1 <- as.vector(unique(as.vector(tHF[tHF$SignHcF1=="*","TaxVar"]))))
(SigVarsHcF2 <- as.vector(unique(as.vector(tHF[tHF$SignHcF2=="*","TaxVar"]))))

# t3 <- t1 %>% 
#   select(TaxVar,VarName,Median) %>% 
#   spread(key = VarName,value = Median, fill = 0) %>% 
#   separate(TaxVar, into = c("Kingdom", "Phylum","Class","Order","Family","Genus"),remove = FALSE,
#            extra = "merge",fill = "left") %>% 
#   mutate(AllVar = ifelse(Genus == "g",paste0(Family,Genus),Genus)) %>% 
#   mutate(AllVar = ifelse(AllVar == "fg", paste0(Order,AllVar), AllVar)) %>% 
#   filter(TaxVar %in% colnames(ydata)[c(280:284)])
# colnames(t3)
# gplots::heatmap.2(as.matrix(t3[,8:19]), scale = "column", col = rev(seriation::bluered(200)), #lhei=lhei, lwid=lwid,lmat =lmat,
#                   trace = "none", density.info = "none")


#Pega os 9 primeiros significativos
SignTax <- as.vector(SigVarsH0$TaxVar)
SignTax <- as.vector(SigVarsH0F1$TaxVar)
SignTax <- as.vector(SigVarsH0F2$TaxVar)
SignTax <- as.vector(SigVarsHcF1$TaxVar)
SignTax <- as.vector(SigVarsHcF2$TaxVar)
SignTax

#Significant Effects within Factor 1
filteredTaxa <- SignTax[3:11]
filteredSign <- tHF %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% #View
  filter(TaxVar %in% filteredTaxa) %>% 
  mutate(SignPosition = ifelse(SignHcF2 == "*", Q3,NA))
summary(filteredSign)

Fact.Res.DF  %>% 
  gather(key = "VarName", value = "Coeff", -TaxVar) %>%
  filter(TaxVar %in% filteredTaxa) %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  left_join(filteredSign,by=c("TaxVar", "Fact1", "Fact2")) %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% #View
  ggplot(aes(x = Fact1, y = Coeff, fill = Fact2)) +
  guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
  stat_summary(fun.y="mean", geom="point", shape=23, size=4, fill="white") +
  stat_smooth(method="lm", formula=y~1, se=FALSE, aes(group = 1)) +
  # stat_summary(fun.y = "mean", color = "red", geom = "line", aes(group = 1))+
  # geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  geom_text(data=filteredSign,aes(x=Fact1,group=Fact2,y=SignPosition,label=SignHcF2),
            size=5, position = position_dodge(width=0.9))+
  # geom_text(aes(label = SignHcF1,x=Fact1, y=Q3),label="*", position = "fill")+
  facet_wrap(~TaxVar,scales = "free_y")+
  theme_bw()

#Significant Effects within Factor 2
filteredTaxa <- SignTax[3:11]
filteredSign <- tHF %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% #View
  filter(TaxVar %in% filteredTaxa) %>% 
  mutate(SignPosition = ifelse(SignHcF1 == "*", Q3,NA))


Fact.Res.DF  %>% 
  gather(key = "VarName", value = "Coeff", -TaxVar) %>%
  filter(TaxVar %in% filteredTaxa) %>% 
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  left_join(filteredSign,by=c("TaxVar", "Fact1", "Fact2")) %>% 
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% #View
  ggplot(aes(x = Fact2, y = Coeff, fill = Fact1)) +
  guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
  stat_summary(fun.y="mean", geom="point", shape=23, size=4, fill="white") +
  stat_smooth(method="lm", formula=y~1, se=FALSE, aes(group = 1)) +
  # stat_summary(fun.y = "mean", color = "red", geom = "line", aes(group = 1))+
  # geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  geom_text(data=filteredSign,aes(x=Fact2,group=Fact1,y=SignPosition,label=SignHcF1),
            size=5, position = position_dodge(width=1))+
  # geom_text(aes(label = SignHcF1,x=Fact1, y=Q3),label="*", position = "fill")+
  facet_wrap(~TaxVar,scales = "free_y")+
  theme_bw()

# Junta todos (faz o inverso)
TaxCode <- taxInfo
Fact.Res.DF  %>%
  filter(TaxVar %in% unique(filteredSign$TaxVar)) %>% #dim()
  gather(key = "VarName", value = "Coeff", -TaxVar) %>% #View()
  # group_by(TaxVar,VarName) %>% 
  # summarise(median = quantiles_95(Coeff)[3]) %>% #View
  separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
  # separate(TaxVar, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), remove = FALSE,
  #          sep = "\\.", fill = "left", extra = "merge") %>% 
  # mutate(AllVar = ifelse(Genus == "g",paste0(Family,Genus),Genus)) %>% 
  # mutate(AllVar = ifelse(AllVar == "fg", paste0(Order,AllVar), AllVar)) %>% #View
  mutate(Fact1 = factor(Fact1, levels = lvF1)) %>% 
  mutate(Fact2 = factor(Fact2, levels = lvF2)) %>% 
  ggplot(aes(x = TaxVar, y = Coeff, fill = Fact2)) +
  guides(fill=F) +
  # geom_boxplot() + 
  stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(limits =fact1)+
  geom_text(data=filteredSign,aes(x=TaxVar,group=Fact2,y=SignPosition,label=SignHcF2),
            size=5, position = position_dodge(width=0.9))+
  facet_wrap(~Fact1, ncol = 2, scales = "free_x")+
  coord_flip()+
  theme_bw()

# Fact.Res.DF.adj  %>%
#   filter(TaxVar %in% SignTax) %>% #dim()
#   gather(key = "VarName", value = "Coeff", -TaxVar) %>% #View()
#   # group_by(TaxVar,VarName) %>% 
#   # summarise(median = quantiles_95(Coeff)[3]) %>% #View
#   separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
#   separate(TaxVar, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), remove = FALSE,
#            sep = "\\.", fill = "left", extra = "merge") %>% 
#   mutate(AllVar = ifelse(Genus == "g",paste0(Family,Genus),Genus)) %>% 
#   mutate(AllVar = ifelse(AllVar == "fg", paste0(Order,AllVar), AllVar)) %>% #View
#   mutate(Fact2 = factor(Fact2, levels = c("7", "14", "21"))) %>% 
#   filter(Fact1 %in% c("RT", "RTA")) %>% 
#   ggplot(aes(x = AllVar, y = Coeff, fill = Fact2)) +
#   guides(fill=F) +
#   # geom_boxplot() + 
#   stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+
#   geom_hline(yintercept = 0, linetype = "dashed")+
#   # scale_x_discrete(limits =fact1)+
#   facet_wrap(~Fact1+Fact2, ncol = 3, scales = "free_x")+
#   coord_flip()+
#   theme_bw()
# 
# # Junta todos (faz o inverso)
# Fact.Res.DF.adj  %>% 
#   gather(key = "VarName", value = "Coeff", -TaxVar) %>% #View()
#   group_by(TaxVar,VarName) %>% 
#   summarise(median = quantiles_95(Coeff)[3]) %>% #View
#   filter(TaxVar %in% SignTax[1]) %>% 
#   separate(VarName, into = c("Fact1", "Fact2"), sep = "_", extra = "merge") %>% 
#   ggplot(aes(x = TaxVar, y = Coeff, fill = Fact2)) +
#   guides(fill=F) +
#   geom_boxplot() + 
#   # stat_summary(fun.data = quantiles_95, geom="boxplot",position = "dodge")+  
#   geom_hline(yintercept = 0, linetype = "dashed")+
#   # scale_x_discrete(limits =fact1)+
#   facet_wrap(~Fact1+Fact2)+
#   coord_flip()+
#   theme_bw()


#Relative abundance

df.abund.predict <- vegan::decostand(selected.model$prediction$ypredMu, method = "total")
heatmap(t(df.abund.predict))
barplot(colSums(df.abund.predict))

library(ggplot2)

data.frame(Fact1=output_all.F$inputs$xdata$Inoculation, #
           Fact2=output_all.F$inputs$xdata$Variety,
           Tax = df.abund.predict[,"Tax100"]) %>%
  ggplot(aes(x = Fact1, y = Tax, fill = Fact2))+
  geom_boxplot()+
  ylab(TaxaName)

data.frame(Fact1=output_all.F$inputs$xdata$Inoculation, #
           Fact2=output_all.F$inputs$xdata$Variety,
           Tax = df.abund.predict[,TaxaName]) %>%
  ggplot(aes(x = Fact1, y = Tax, fill = Fact2))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean")+
  ylab(TaxaName)
