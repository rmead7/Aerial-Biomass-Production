###Biomass Production in the Cape Fear Estuary of North Carolina.###
##Part I:##
#Loading Dataframe
df <- read.table("/Users/robbiemead/Documents/Graduate School/Courses/Fall 2021/MATH 564/Project /Project/LINTHALL.txt",header = TRUE)
df <- subset(df,select = -c(Obs,Loc,Type))
View(head(df))

#OLSR
model <- lm(BIO~H2S+SAL+Eh7+pH+BUF+P+K+Ca+Mg+Na+Mn+Zn+Cu+NH4,data = df)
coef <- as.data.frame(model$coefficients)
coef
varnum <- as.data.frame(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15"))
coef <- cbind(varnum,coef)
coef
#Collinearity Diagnostics
#i.)standardized residuals versus fitted values
library('ggplot2')
library("ggthemes")
ggplot(model, aes(x = model$fitted.values,y = model$residuals)) + 
  geom_point()+theme_minimal()+labs(x = "Fitted Values",y = "Residuals",title = "Predicted v. Residuals")

plot(model$fitted.values,model$residuals,xlab = "Predicted",ylab = "Residuals",main = "Predicted v. Residuals")
abline(h = 0,col = "red",lwd = 3,lty = 3)

#ii.)pairwise correlation coefficient matrix
install.packages('kable')
library('kable')
pcor <- cor(df)[-1,-1]
pcorr <- round(pcor,3)
View(pcorr)
pcorr <- print(pcor,digits = 4)
kable(pcorr)
pairs(df,main = "Pairwise Correlation Matrix")

#pH and BUF have a high correlation of -0.946, pH and Ca have a high correlation of 0.8779, pH and Zn have a high correlation of -0.722
#pH and NH4 have a high correlation of -0.745, BUF and Ca have a high correlation of -0.791, BUF and Zn have a high correlation of 07.14,
#BUF and Zn have a high correlation of 0.849, K and Mg have a high correlation of 0.862, K and Na have a high correlation of 0.792, 
#K and Cu have a high correlation of 0.693, Ca and Zn have a high correlation of -0.699, Mg and Na have a high correlation of 0.899,
#Mg and Cu have a high correlation of 0.712, Zn and NH4 have a high correlation of 0.720.

#iii.) VIF  
install.packages("regclass")
library('regclass')

VIF <- as.data.frame(VIF(model))
colnames(VIF)[1] <- "VIF"
#VIF shows that the predictor variables pH, BUF, Ca, Mg, Na, Zn are heavily affected by collinearity, because
#the VIF values exceed 10. The predictor variables H2S, SAL, K,Mn, CU and NH4 show that they are also effected
#by collinearity, since there VIF is between the interval of 3 and 10. The predictor values of Eh7 and P are showing
#promise of no effect of collinearity. 

##Part II:##
#Principle Component Regression
#i.)standardizing the dataframe
bio <- (df$BIO - mean(df$BIO))/sd(df$BIO)
h2s <- (df$H2S - mean(df$H2S))/sd(df$H2S)
sal <- (df$SAL - mean(df$SAL))/sd(df$SAL)
eh7 <- (df$Eh7 - mean(df$Eh7))/sd(df$Eh7)
ph <- (df$pH - mean(df$pH))/sd(df$pH)
buf <- (df$BUF - mean(df$BUF))/sd(df$BUF)
p <- (df$P - mean(df$P))/sd(df$P)
k <- (df$K - mean(df$K))/sd(df$K)
ca <- (df$Ca - mean(df$Ca))/sd(df$Ca)
mg <- (df$Mg - mean(df$Mg))/sd(df$Mg)
na <- (df$Na - mean(df$Na))/sd(df$Na)
mn <- (df$Mn - mean(df$Mn))/sd(df$Mn)
zn <- (df$Zn - mean(df$Zn))/sd(df$Zn)
cu <- (df$Cu - mean(df$Cu))/sd(df$Cu)
nh4 <- (df$NH4 - mean(df$NH4))/sd(df$NH4)
dfs <- data.frame(BIO = bio, H2S = h2s, SAL = sal, EH7 = eh7, pH = ph, BUF = buf, P = p, K = k,
                  Ca = ca, Mg = mg, Na = na, Mn = mn, Zn = zn, Cu = cu, NH4 = nh4)

#ii.)pairwise correlation coefficient matrix of standardized data frame
pscor <- cor(dfs)
View(pscor)

#iii.)eigenvalues and eigenvectors of standardized correlation coefficient matrix
eigenvectors_std <- eigen(pscor)
eigenvectors_std <- eigenvectors_std$vectors
eigenvectors_std
eigenvalues_std <- eigen(pscor)
eigenvalues_std <- eigenvalues_std$values
eigenvalues_std

#iv.)construct principle component data frame with BIO from the original data frame. 
BIOX <- dfs[,-1]
BIOX.pca <- prcomp(BIOX,center = TRUE,scale. = TRUE)
BIOX.pca$rotation
biomass.pcr.df <- cbind(df[,1],data.frame(BIOX.pca$x))
colnames(biomass.pcr.df)[1] <- "BIO"
head.biomass <- head(biomass.pcr.df)
View(biomass.pcr.df)

#v.)Correlation between PC's and BIO
cor(biomass.pcr.df)[,1]
#PC6, PC10, PC12, PC13, PC14 have a very little correlation with BIO since there values are negligable. 

#vi.)full regression model
model.pcr <- lm(BIO~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14,data = biomass.pcr.df)
sum <- summary(model.pcr)
#There are a significant amount of PC's that are considered not significant in the full regression model - the following PC's that can be excluded from
#the PC model are the following: PC4, PC5, PC6, PC10,PC11,PC12,PC13,PC14. The reduced PC regression model can be explained below.


#viii.)computing the beta values from the principle component regression reduced model. 
betas <- BIOX.pca$rotation %*% model.pcr$coefficients[-1] 
betas <- as.data.frame(betas)
colnames(betas)[1] <- "Beta Values"
betas

rmcoef <- model.pcr.reduced$coefficients
rmcoef["PC4"]=0
rmcoef <- c(1000.8,216.1552,78.8904,133.1758,0,0,0,0,194.4449,330.9289,0,0,0,0,0)
betasr <- BIOX.pca$rotation %*%as.vector(rmcoef)[-1]
betasr

##Part III:##
#Loading Dataframe
df_3 <- read.table("/Users/robbiemead/Documents/Graduate School/Courses/Fall 2021/MATH 564/Project /Project/LINTH-5.txt",header = TRUE)
df_3 <- subset(df_3,select = -c(Obs,Loc,Type))
View(df_3)

#Stepwise Regression
#i.)step 1: BIO~SAL, BIO~pH, BIO~K, BIO~Na, BIO~Zn
anova(lm(BIO~SAL,df_3))
summary(lm(BIO~SAL,df_3))
summary(lm(BIO~pH,df_3))
summary(lm(BIO~K,df_3))
summary(lm(BIO~Na,df_3))
summary(lm(BIO~Zn,df_3))
#Since pH has the smallest p-value, that predictor will remain. 

#ii.) step 2: BIO~pH + SAL, BIO~pH + K, BIO~pH + Na, BIO~pH + Zn
summary(lm(BIO~pH+SAL,df_3))
summary(lm(BIO~pH+K,df_3))
summary(lm(BIO~pH+Na,df_3))
summary(lm(BIO~pH+Zn,df_3))
#Since the p-values remained accepted at the alpha value of alpha = 0.15 for both pH and Na,
#we can accept x7 into the regression model. The p-values for the regression of Na were the smallest,
#while maintaining the integrity of pH. 

#iii.)step 3: y~x4 + x10 + x2, y~x4 + x10 + x7, y~x4 + x10 + x12
summary(lm(BIO~pH+Na+SAL,df_3))
summary(lm(BIO~pH+Na+K,df_3))
summary(lm(BIO~pH+Na+Zn,df_3))
#Since the introduction of any of the other predictors are not significant, then, the model we have currently
#is the model that is the best model through stepwise regression. 

#Verifying Stepwise Regression
install.packages('leaps')
library('leaps')
library('MASS')

full.model <- lm(BIO~.,df_3)
stepwise.model <- stepAIC(full.model,direction = "both")
summary(stepwise.model)

#The stepwise regression shows that the best model is when alpha is 0.15, there are four predictors, with 
#SAL, pH, K, and Zn because there alpha values are below 0.15. So the final model that shows a significant result 
# would be classified as BIO =  1505.4882 - 35.9433SAL + 293.8611pH - 0.4388K - 23.4519Zn. 

#Subset Selection
install.packages('caret')
library('caret')

#Creating the subset model
step.model <- regsubsets(BIO~.,data = df_3,nbest = 2,method = "exhaustive")
subset.summary<- summary(step.model)
subset.summary
subset.summary.outmat <- as.matrix(subset.summary$outmat)

#Identifying vital statistics to make conclusion
Cp <- as.vector(subset.summary$cp)
VIF <- as.vector((1)/(1-subset.summary$rsq))
best.subset <- cbind(subset.summary.outmat,Cp,VIF)
View(best.subset)

#Though assessing the best sub set for each of the models that have 2 predictors, the best model there is 
# when pH and Na are in the model because the Cp value, Cooks Distance, is smaller than that of the other model
#pH and K because the Cp value is larger, bt only by small margins. Interestingly, the VIF values are switched where
#the VIF of the second best 2 predictor model is smaller than that of the best 2 predictor model, but again, only
#by small margins. 
