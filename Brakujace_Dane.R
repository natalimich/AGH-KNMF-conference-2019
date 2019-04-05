#www.ubs.com/polandcareers

library(mice)
library(VIM)
library(BaylorEdPsych)
library(plyr)
library(dplyr)
library(ggplot2)

estimationDataNA <- read.csv("FinancialData2.csv")
CriteriaSummary <- read.csv("CriteriaSummary.csv")

ChampionModelCriteria <- c('g1_15_AQ','g1_25_AQ','g2_38_CL','g2_37_CL','g4_61_LF','g4_68_LF','g5_86_PE','g5_99_PE')

VisualizationData<- estimationDataNA[,c(ChampionModelCriteria)]

# abbreviations
dataDescription <- CriteriaSummary[,c("Code","Name")]

variableNames <- data.frame(merge(dataDescription, data.frame(code = names(VisualizationData)), by.x="Code", by.y="code", all.y = FALSE))
variableNames_names <- as.factor(variableNames[,2])
variableNames_codes <- as.factor(variableNames[,"Code"])

abbreviationsFun<- function(varNames){
  abbreviations <- sapply(strsplit(x=as.character(varNames),split="[\\+\\\\$\\ \" \"-]"),function (x) paste(sapply(x, function (y) if (toupper(y) == y) {y} else {if (y %in% c("and", "of", "the", "to")) {substr(y,1,1)} else {toupper(substr(y,1,1))}}),collapse="")) 
  
  for (i in 1:length(abbreviations)){               
    abbreviations[i]<-gsub("\\n","",abbreviations[i],perl=TRUE)
    #abbreviations[i]<-gsub("\\+","",abbreviations[i],perl=TRUE)
    abbreviations[i]<-gsub("\\((?!C|N|L|%)", "",abbreviations[i],perl=TRUE)
    abbreviations[i]<-gsub("(?<!C|N|L|%|<)\\)","",abbreviations[i],perl=TRUE)
    parted<-unlist(strsplit(abbreviations[i], "(?=[[:punct:]])", perl=TRUE))
    unclosed<-length(grep("\\(",parted))-length(grep("\\)",parted))
    unopened<-length(grep("\\)",parted))-length(grep("\\(",parted))
    if(unclosed>0){
      parted<-parted[-grep("\\(",parted)[1:unclosed]]
    }
    if(unopened>0){
      parted<-parted[-grep("\\)",parted)[1:unopened]]
    }
    abbreviations[i]<-paste(parted,collapse="")
  }
  return(abbreviations) 
}

varAbbreviations <- abbreviationsFun(variableNames_names)
#varAbbreviations <- cbind(as.charactehead(r(variableNames_codes), varAbbreviations)
#Pattern-plot my code 

PatternPlot <- function(dataToVis, abbreviations){
  dane <- data.matrix(dataToVis)
  r <- 1 * is.na(dane) 
  n.var <- ncol(dane) 
  n.obs <- nrow(dane) 
  mdp <- (r %*% (2^((1:n.var - 1)))) + 1  # sing a number of pattern to split the data into groups 
  r.matrix <- cbind(r,mdp) #assing the number of pattern to the pattern 
  r.o.matrix <- as.data.frame(r.matrix[order(mdp),]) #ordering by pattern number 
  colnames(r.o.matrix) <- c(abbreviations, "mdp") 
  freq <- data.frame(count(r.o.matrix, "mdp")) 
  r.o.matrix <- r.o.matrix[order(-colSums(r.o.matrix))]
  miss.per <- paste0(round((colSums(r.o.matrix[,!(names(r.o.matrix) %in% c("mdp","freq","n")) ])*100/n.obs),2),"%")
  new.label <- NULL
  p = 1
  for( i in 1:length(abbreviations)){
    new.label[p]<- paste0(colnames(r.o.matrix[,!(names(r.o.matrix) %in% c("mdp","freq","n")) ])[i]," (",miss.per[i],") " )
    p=p+1
  }
  pattern.matrix2 <- merge(r.o.matrix, freq, by = "mdp", all.x=TRUE)
  pattern.matrix2[,"freq"]<- round((pattern.matrix2[,"freq"]*100)/n.obs,2)
  pattern.matrix2 <- unique(pattern.matrix2)
  pattern.matrix2 <- (pattern.matrix2[order(-pattern.matrix2[,"freq"]),])
  freq.per <- paste0(pattern.matrix2[,"freq"], "%")
  pattern.to.plot <- pattern.matrix2[,!(names(pattern.matrix2) %in% c("mdp","freq","n")) ]
  m <- as.matrix(pattern.to.plot)
  t.mat <- t(m)
  t.mat.plot=t.mat
  colnames(t.mat.plot) <- 1:ncol(t.mat.plot)
  x11=melt(t.mat.plot)
  names(x11)=c("x","y","Variable")
  x11$Variable=factor(x11$Variable)
  
  levels(x11$Variable)=c("observed","missing")
  #intresting
  print(qplot(x, y, fill=Variable, data=x11, geom='tile')+ scale_fill_manual(values = c("red","blue"))+ scale_y_continuous("No. Pattern",breaks=seq(1.5,length(freq$mdp)+0.5,1),label=seq(1,length(freq$mdp),1),sec.axis=dup_axis(name="Frequency", labels=c(freq.per)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10,color="black",margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),axis.text.y = element_text(angle = 0, hjust = 1, vjust = 1,size=8,color="black" ))+theme(axis.title.y = element_text(size = rel(1), angle = 90))+ theme(axis.title.x = element_text(size = rel(1), angle = 0,vjust=3))+coord_cartesian(expand = FALSE) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + labs(x="Model Variables")+geom_vline(xintercept = 0.5:(length(varAbbreviations)+1))+ geom_hline(yintercept = 0.5:(length(freq$mdp)+1)) + scale_x_discrete(labels=new.label))
}

PatternPlot(VisualizationData, varAbbreviations)

#get colnames to Visualization data 
colnames(VisualizationData) <- varAbbreviations
aggr(VisualizationData, comined=TRUE)

#MATRIXPLOTS

matrixplot(VisualizationData, sortby = varAbbreviations[2])


#SCATTERPLOTS

scattmatrixMiss(log10(VisualizationData), highlight = c(varAbbreviations[2]),col=c("grey","black"), pch=16,interactive=FALSE)


#HISTOGRAMS
varAbbreviationsHIST <- colnames(VisualizationData)

getPairs <- function(positon.horizontal){
  horizontal <- (varAbbreviationsHIST[positon.horizontal])
  horizontal <- rep(horizontal, (length(varAbbreviationsHIST)-1))
  vertical <- (varAbbreviationsHIST[-(positon.horizontal)])
  a <- cbind(horizontal, vertical)
  colnames(a) <- NULL 
  return(a)
}
  par(mfrow=c(round((length(varAbbreviationsHIST)-1)/2),2))
  
  for( i in 1:(length(varAbbreviationsHIST)-1)){
  histMiss(VisualizationData[,c(getPairs(2)[i,])],interactive=FALSE, col=c("white","black"))
}


####Imputation####
methods(mice) 

dataNAomit <- na.omit(estimationDataNA)
modelFit <- glm(deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_PE, data=dataNAomit, family=binomial())
summary(modelFit)


#mean
meanData <- mice(estimationDataNA, method="mean")
#print(xyplot(meanData, deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_PE, cex=1))
densityplot(meanData)
#xyplot(meanData, deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_PE,pch=18,cex=1)
stripplot(meanData, pch = 20, cex = 1.2)
modelFit1 <- with(meanData,glm(deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_PE, family = binomial(link='logit')))
pooled1 <- pool(modelFit1)
summary(pooled1)


#pmm
pmm1Data <- mice(estimationDataNA,m=5,maxit=50,meth='pmm',seed=500)
#xyplot(pmm1Data, deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_PE, cex=1)
densityplot(pmm1Data)
stripplot(pmm1Data, pch = 20, cex = 1.2)
plot(pmm1Data)
pmm1Data$pred
modelFit2 <- with(pmm1Data,glm(deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_PE,  family=binomial(link='logit')))
pooled2 <- pool(modelFit2)
summary(pooled2)[,c("est","fmi","lambda")]



#pmm2
pmm2Data <- mice(estimationDataNA,m=50,maxit=50,meth='pmm',seed=500)
#xyplot(pmm2Data, deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_P, cex=1)
densityplot(pmm2Data)
stripplot(pmm2Data, pch = 20, cex = 1.2)
plot(pmm2Data)
modelFit3 <- with(pmm2Data,glm(deflag~ g1_15_AQ+g1_25_AQ+g2_38_CL+g4_61_LF+g4_68_LF+g5_86_PE+g5_99_PE,family=binomial(link='logit')))
pooled3 <- pool(modelFit3)
summary(pooled3)






