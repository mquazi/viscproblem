Infineum PROBLEM

install.packages("stargazer")
install.packages("texreg")
install.packages("leaps")
install.packages("car")
install.packages("MASS")
install.packages("lmerTest")
install.packages("lmtest")
install.packages("lme4")
library(lme4)
library(lmtest)
library(leaps)
library(car)
library(leaps)
library('ggplot2')
library(stargazer)
library(texreg)
library(MASS)
library(car)
library(dae)
library(lmerTest)

data<-read.csv(file.choose(),header=T)
summary(data)
attach(data)
table(data$DI)
table(data$VM)

tapply(Blend.KV40,DI,max)
which.max(Blend.KV40)
summary(Blend.KV40,DI)
table(data$VM)
names(data)
nrow(data)
attach(data)

boxplot(Blend.KV40~DI,ylab="Blend.KV40", xlab="DI", main = 
                "Boxplot BlendKV40 vs DI")
boxplot(Blend.KV40~VM,ylab="Blend.KV40", xlab="VM", main = 
                "Boxplot BlendKV40 vs VM") 



##Continuous variables only
contdata<-cbind(Blend.KV40,DI,VM,BS.Density,BS.KV100,BS.KV40, BS.total)
head(contdata)
pairs(contdata)
cor(contdata)



b) Model selection
y<-Blend.KV40
y
x1<-DI
x2<-VM
x3<-BS.Density
x4<-BS.KV100
x5<-BS.KV40
x6<-BS.total
plot(x1*x2,y)
cor(x1*x2,y)
0.55
cor(x1*x5,y)
0.53
cor(x2*x3,y)
0.67
cor(x2*x4,y)
0.90
cor(x2*x5,y)
0.94
cor(x2*x6,y)
0.63
cor(x3*x4,y)
0.71
cor(x3*x5,y)
0.71
cor(x4*x5,y)
0.71
cor(x5*x6,y)
0.61

# Most promising interactions for later consideration... 
x2x4<-x2*x5
x2x5<-x2*x5



maindata<-cbind(y,x1,x2,x3,x4,x5,x6,
                x2x4,x2x5)
head(maindata)
pairs(maindata)
cor(maindata)
cordata<-cbind(x2,x4,x5,x2x4,x2x5)

cor(cordata)


##deleted x2 interactions because of high cor

##Backward elimination

upper<-formula(~x1+x2+x3+x4+x5+x6)
lower<-formula(~1)

mod1full<-lm(y~x1+x2+x3+x4+x5+x6)
mod1step<-step(mod1full,scope=list(lower=lower,upper=upper),direction="backward")
summary(mod1step)
backmodel<-lm(y~x1+x2+x3+x5)
vif(backmodel)

#### Forward
upper<-formula(~x1+x2+x3+x4+x5+x6)
lower<-formula(~1)
class(maindata)
data<-as.data.frame(maindata)
mod1null<-lm(y~1,data=data)
mod1forward<-step(mod1null,scope=list(lower=lower,upper=upper),direction="forward")
summary(mod1forward)

## Stepwise
upper<-formula(~x1+x2+x3+x4+x5+x6)
lower<-formula(~1)

mod1<-lm(y~1,data=data)
mod1stepwise<-step(mod1full,scope=list(lower=lower,upper=upper),direction="both")
summary(mod1stepwise)


### BEST SUBSET

library(leaps)
## Cp
X<-data.frame(cbind(x1,x2,x3,x4,x5,x6))
X<-as.matrix(X)
y
cpsubset<-leaps(X,y,method="Cp",nbest=2)
cpsubset$which[order(cpsubset$Cp)[1:5],]
cpsubset$Cp[order(cpsubset$Cp)[1:5]]
cpsubset
cpsubset$Cp

##r2, adjr2, cp, bic
maindata<-data.frame(maindata)
leaps=regsubsets(y~x1+x2+x3+x4+x5+x6,data=maindata, nbest=2)
plot(leaps, scale="adjr2")
plot(leaps, scale="Cp")
plot(leaps, scale="bic")
plot(leaps, scale="r2")
summary(leaps)
leaps

###Final Model selected is
finalmodel<-lm(y~x1+x2+x3+x5)
summary(finalmodel)



##Diagnostics

plot(finalmodel)
plot(fitted.values(finalmodel),resid,main="Residuals vs Fitted values") ##major problem
resid<-residuals(finalmodel)
sum(resid)

#Non const var
bptest(y~x1+x2+x3+x5)
p>0.05, p=0.2545
const error var form the test, but the test is failing here... 

#Normality
par(mfrow=c(1,1))
hist(resid)
boxplot(resid,main="Boxplot residuals")
shapiro.test(resid)
p<alpha, p=1.214e-05
not normal... 

## Outliers
X outliers
leverage<-hatvalues(finalmodel)
sum(leverage)
xoutliers<-which(leverage > 3*5/86)
xoutliers
plot(leverage)

Y outliers
rstud<-rstudent(finalmodel)
outlierTest(finalmodel)
p=0.067602, p>alpha
No outlier... 

rstuden<-abs(rstud)
youtliers<-which(rstuden >= qt(1-0.01/(2*86),80))
youtliers
plot(cooks.distance(finalmodel),pch=23,main="Cook's distance")
No outlier... 


##
vif(finalmodel)

##transformation

library(MASS)
boxcox(y~x1+x2+x3+x5,lambda=seq(-2,2,length=10))
gmean<-exp(mean(log(y)))
sse<-c()
lambda<-c()
i<-1
for (lam in seq(-2,2,0.1)){
        if (lam !=0){
                k1<-(y^lam - 1)/(lam*gmean^(lam-1))
        } else {
                k1<-log(y)*gmean
        }
        test<-anova(lm(k1~x1+x2+x3+x5))
        sse[i]<-test['Residuals','Sum Sq']
        lambda[i]<-lam
        i<-i+1
}
cbind(lambda,sse)

abline(0,mean(sse),lty=2)
abline(mean(sse),0,lty=2)

ynew<-y^(0.5)
y
ynew
transformedmodel<-lm(ynew~x1+x2+x3+x5)
summary(transformedmodel)

##Diagnostics for transformed

plot(transformedmodel)
transresid<-residuals(transformedmodel)

plot(fitted.values(transformedmodel),transresid,main="Transformed Model Residuals vs Fitted values") ##major problem

#Non const var
bptest(ynew~x1+x2+x3+x5)
p>0.05, p=0.1442
const error var but curvature is much better than before, eventhough I feel the test is 
# is failing here 

#Normality
par(mfrow=c(1,1))
hist(transresid)
boxplot(transresid,main="Boxplot Transformed model residuals")
shapiro.test(transresid)
p<alpha, p=0.0002626 ##Not normal but major improvement

## Outliers
X
leverage<-hatvalues(transformedmodel)
sum(leverage)
plot(leverage)
xoutliers<-which(leverage > 3*5/86)
xoutliers
plot(leverage)

Y
rstud<-rstudent(transformedmodel)
outlierTest(transformedmodel)
rstuden<-abs(rstud)
youtliers<-which(rstuden >= qt(1-0.01/(2*86),80))
youtliers
plot(cooks.distance(transformedmodel),pch=23,main="Cook's distance")


Dffits
dffits(transformedmodel)
2*(sqrt(4/86))
maindata[which(dffits(transformedmodel)>0.4313311),]
maindata[which(dffits(finalmodel)>0.4313311),]
plot(dffits(transformedmodel),pch=2,ylab="DFFITS", xlab="Obs no.", main="DFFITS")
plot(dffits(finalmodel),pch=2,ylab="DFFITS", xlab="Obs no.", main="DFFITS")


##
vif(transformedmodel)
#####added variable plots
plot(transresid,x2x4)
plot(transresid,x2x5)


summary(transformedmodel)
model<-lm(ynew~x1+x2+x3+x5)

summary(model)

plot(model)



Tables

model2realstuff<-lm(ynew~x1+x2+x3+x5)
summary(model2realstuff)
stargazer(anova(model),summary=F, title="ANOVA Table for Final Model")
stargazer((model2realstuff),summary=F)
stargazer(vif(model,finalmodel),title="Variance inflation Factors Final Model")
update.packages("stargazer")
library(stargazer)

