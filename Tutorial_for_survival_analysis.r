# Install require packages
install.packages("ggpubr")
install.packages("ggplot2")
install.packages("survival")
install.packages("prodlim")
# Load require packages 
library(ggplot2)
library(ggpubr)
library(survival)
library(prodlim)

######1. Draw boxplot
# Set the environment where the files are located
setwd("C:/Users/nttg8/Desktop/BRCA/BRCA/Figure/Boxplot")

data=read.delim("tenfile.txt",header=TRUE)
p2 <- ggplot(data, aes(x=Genes, y=EX, fill=BCLM),ylab="log2 gene expression at mRNA level") +geom_boxplot() +facet_wrap(~Genes, scale="free")
p2+ylab("Log2 gene expression at mRNA level")+stat_compare_means(label =  "p.signif", label.x = 1.5)


######2. Draw survival files
# Set environment
setwd("C:/Users/Mai Nhu Nguyet/Downloads/test")
V <- read.delim("gse2034 sur.txt")

# Create the output folder with pdf file inside
pdf(file.path("output", "KMP GSE2034 52 probs.pdf"))
# Draw survival curve using hist function on prodlim package
km1 <- prodlim(Hist(mBM,BM)~Risk,data=V)
# Set up the paremeter for the plot
par(mar=c(4,1,1,1), mgp=c(3,1,0), cex.lab=1.4, font.lab=1, font.axis=2,  ps=16, family="Helvetica")
plot(km1, marktime=TRUE, marktime.cex=0.7, col=c("red", "blue"),
     atrisk.labels=paste(c("High risk   ","Low risk    ")),
     axis2.las=2, ylim=c(0,1), axis2.at=seq(0,1,0.20), axis2.labels=c(0,20,40,60,80,100),
     atrisk.title="", lwd=1.5, #line width for all curves
     atrisk.at=seq(0,200,20), #time grid for no. At risk
     ylab="Probability of Bone metastasis (%)", xlab="Metastasis (months)",  # label for x-axis
     axis1.at=seq(0,200,20), # time grid for x-axis
     legend=FALSE,
     atrisk.cex=1,
     atrisk.adj=1.3,
     atrisk.line=c(0.5,2),
     logrank=FALSE, background=FALSE, confint.citype=FALSE) # show log-rank p-value if change to TRUE
legend (90, 0.15, c("High risk", "Low risk"), lty = c(1,1), col = c ("red","blue"), bty = "n", y.intersp = 1.7)
y=1- pchisq(survdiff(Surv(mBM,BM) ~Risk, data=V)$chisq, 1)
text(20,0.015,paste("p=",signif(y,3)), font =3) 
text(-3,-.2,"No. at risk",xpd=NA, family="Helvetica", font = 2, font = 2) 
dev.off()

#######3. Export the test  of chisquare, fisher exact and  univariate and multivarite cox regression analysis test
# set enviroment
setwd("Dien cai duong dan vao day nhe !")

#Removing rows with missing values in the file V, the file shoud be named after underscore ___ instead of space " "
V <- read.delim("asmt exp 06092019.txt")
attach(V)
V0 <- read.delim("asmt exp 06092019.txt")
attach(V0)

V <- na.omit(V0)
d.complete <- na.omit(d)

# Check for the chisquare and fisher exact test between our model and the clinical information
summary(V)
table(age51,Risk)
# To transform data to percent
t<-prop.table(table(age51,Risk))
round(t*100,1)
chisq.test(table(age51,Risk),simulate.p.value = TRUE)
fisher.test(table(age51,Risk))

# Coxregression analysis
coxph(formula=Surv(mBone,Bone) ~ BCR, data=V)
summary(coxph(formula=Surv(mBone,Bone) ~ BCR, data=V))
summary(coxph(formula=Surv(mBone,Bone) ~ BCR + factor(ESR1) +  factor(PR) + factor(HER2), data=V))
summary(coxph(formula=Surv(mBM,BM) ~ ESR1, data=V))
summary(coxph(formula=Surv(mBone,Bone) ~ BCR + ESR1 +  PR + HER2+Age, data=V))
