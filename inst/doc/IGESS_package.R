### R code from vignette source 'IGESS_package.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: IGESS-prelim
###################################################
library("IGESS")


###################################################
### code chunk number 3: IGESSExample-prelim
###################################################

data(IGESSDB)
X <- Model$X
y <- Model$y
AIPVal <- Model$AIPVal
dim(AIPVal)
dim(X)


###################################################
### code chunk number 4: IGESS-noZA-show (eval = FALSE)
###################################################
## fit <- IGESS(X, y)


###################################################
### code chunk number 5: IGESS-noZA-show (eval = FALSE)
###################################################
## fit <- iGess(X, y, SS = AIPVal)


###################################################
### code chunk number 6: IGESS-noZA-show (eval = FALSE)
###################################################
## fit <- iGessPlink(genoplinkfile)


###################################################
### code chunk number 7: IGESS-noZA-show (eval = FALSE)
###################################################
## fit <- iGessPlink(genoplinkfile, summaryfile, configfile)


