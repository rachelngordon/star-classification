
# read in the data
dat <- read.csv("Skyserver_SQL2_27_2018 6_51_39 PM.csv", header = TRUE)

# separate out a test set for evaluating accuracy
set.seed(1234)
ind <- sample(1:10000, 3000)
test <- dat[ind,]
SDSS <- dat[-ind,]

# summary of the data
summary(SDSS)

# visualize the frequency of each class within the data
class_freq <- table(SDSS$class)
barplot(class_freq, xlab = "Frequency", ylab = "Class", main = "Class Frequency", 
        col = "cornflowerblue", border = NA, horiz = TRUE)
# most observations were galaxies and stars, not very many quasars

# convert the observed classes to numerical levels
SDSS[SDSS=="STAR"] <- 0
SDSS[SDSS=="GALAXY"] <- 1
SDSS[SDSS=="QSO"] <- 2
SDSS$class <- as.integer(SDSS$class)

# check for missing data
is.na(SDSS)

# correlation matrix for all the variables in the data set
cor(SDSS)

# fit a baseline category logit model with all possible explanatory variables
# check for significance of predictors
library(VGAM)
fit <- vglm(class ~ ra + dec + u + g + r + i + z + run + camcol + field + redshift + plate 
            + mjd + fiberid, family = multinomial, data = SDSS)
summary(fit)

# check for significance when used as sole predictors
# remove variables with p-values > 0.2
fit.ra <- vglm(class ~ ra, family = multinomial, data = SDSS)
summary(fit.ra)
fit.dec <- vglm(class ~ dec, family = multinomial, data = SDSS)
summary(fit.dec)
fit.i <- vglm(class ~ i, family = multinomial, data = SDSS)
summary(fit.i)
fit.r <- vglm(class ~ r, family = multinomial, data = SDSS)
summary(fit.r)
fit.run <- vglm(class ~ run, family = multinomial, data = SDSS)
summary(fit.run)
fit.camcol <- vglm(class ~ camcol, family = multinomial, data = SDSS)
summary(fit.camcol)
fit.field <- vglm(class ~ field, family = multinomial, data = SDSS)
summary(fit.field)
fit.plate <- vglm(class ~ plate, family = multinomial, data = SDSS)
summary(fit.plate)
fit.mjd <- vglm(class ~ mjd, family = multinomial, data = SDSS)
summary(fit.mjd)
fit.fiberid <- vglm(class ~ fiberid, family = multinomial, data = SDSS)
summary(fit.fiberid)
# not significant: camcol, field

# initial main effects model
fit.me <- vglm(class ~ ra + dec + u + g + r + i + z + run + redshift + plate + mjd + fiberid, 
               family = multinomial, data = SDSS)
summary(fit.me)

# conduct backward elimination
AIC(fit.me)
fit1 <- vglm(class ~ ra + u + g + r + i + z + run + redshift + plate + mjd + fiberid, 
             family = multinomial, data = SDSS)
AIC(fit1)
fit2 <- vglm(class ~ ra + u + g + r + i + z + redshift + plate + mjd + fiberid, 
             family = multinomial, data = SDSS)
AIC(fit2)
fit3 <- vglm(class ~ ra + u + g + r + z + run + redshift + plate + mjd + fiberid, 
             family = multinomial, data = SDSS)
AIC(fit3)
fit4 <- vglm(class ~ ra + u + g + r + z + run + redshift + plate + fiberid, 
             family = multinomial, data = SDSS)
AIC(fit4)
fit5 <- vglm(class ~ ra + u + g + z + run + redshift + plate + fiberid, 
             family = multinomial, data = SDSS)
AIC(fit5)
fit6 <- vglm(class ~ ra + u + g + z + run + redshift + plate, 
             family = multinomial, data = SDSS)
AIC(fit6)
# AIC is larger by almost 20
fit7 <- vglm(class ~ ra + u + g + z + redshift + plate, 
             family = multinomial, data = SDSS)
AIC(fit7)

# add additional significant variables back into the model
fit.dec2 <- vglm(class ~ dec + ra + u + g + z + run + redshift + plate, 
                 family = multinomial, data = SDSS)
summary(fit.dec2)
fit.r2 <- vglm(class ~ r + ra + u + g + z + run + redshift + plate, 
                 family = multinomial, data = SDSS)
summary(fit.r2)
fit.i2 <- vglm(class ~ i + ra + u + g + z + run + redshift + plate, 
               family = multinomial, data = SDSS)
summary(fit.i2)
fit.camcol2 <- vglm(class ~ camcol + ra + u + g + z + run + redshift + plate, 
               family = multinomial, data = SDSS)
summary(fit.camcol2)
fit.field2 <- vglm(class ~ field + ra + u + g + z + run + redshift + plate, 
               family = multinomial, data = SDSS)
summary(fit.field2)
fit.mjd2 <- vglm(class ~ mjd + ra + u + g + z + run + redshift + plate, 
                   family = multinomial, data = SDSS)
summary(fit.mjd2)
fit.fiberid2 <- vglm(class ~ fiberid + ra + u + g + z + run + redshift + plate, 
                   family = multinomial, data = SDSS)
summary(fit.fiberid2)
# none significant

# check for plausible interactions
fit.rau <- vglm(class ~ ra + u + g + z + run + redshift + plate + ra*u, 
                family = multinomial, data = SDSS)
summary(fit.rau)
fit.rag <- vglm(class ~ ra + u + g + z + run + redshift + plate + ra*g, 
                family = multinomial, data = SDSS)
summary(fit.rag)
fit.raz <- vglm(class ~ ra + u + g + z + run + redshift + plate + ra*z, 
                family = multinomial, data = SDSS)
summary(fit.raz)
fit.ug <- vglm(class ~ ra + u + g + z + run + redshift + plate + u*g, 
                family = multinomial, data = SDSS)
summary(fit.ug)
fit.uz <- vglm(class ~ ra + u + g + z + run + redshift + plate + u*z, 
                          family = multinomial, data = SDSS)
summary(fit.uz)
fit.gz <- vglm(class ~ ra + u + g + z + run + redshift + plate + g*z, 
                family = multinomial, data = SDSS)
summary(fit.gz)
fit.ugz <- vglm(class ~ ra + u + g + z + run + redshift + plate + u*g*z, 
                family = multinomial, data = SDSS)
summary(fit.ugz)
fit.rrs <- vglm(class ~ ra + u + g + z + run + redshift + plate + run*redshift, 
                family = multinomial, data = SDSS)
summary(fit.rrs)
fit.rp <- vglm(class ~ ra + u + g + z + run + redshift + plate + run*plate, 
               family = multinomial, data = SDSS)
summary(fit.rp)
# very significant
fit.rsp <- vglm(class ~ ra + u + g + z + run + redshift + plate + redshift*plate, 
                family = multinomial, data = SDSS)
summary(fit.rsp)

# final model
fit.final <- vglm(class ~ ra + u + g + z + run + redshift + plate + redshift*plate, 
                  family = multinomial, data = SDSS)
summary(fit.final)

# check the significance using the likelihood ratio test
fit0 <- vglm(class ~ 1, family = multinomial, data = SDSS)
lrtest(fit6, fit0)
lrtest(fit.final, fit6)
# model is very significant without the interaction term
# interaction term is also very significant on its own


# binomial logistic regression model for predicting if an observation is a star or not
SDSS[SDSS==2] <- 1
binom.fit <- glm(class ~ ra + dec + u + g + r + i + z + run + camcol + field + redshift 
            + plate + mjd + fiberid, family = binomial, data = SDSS)
# check for significance
summary(binom.fit)

# check for significance when used as sole predictors
# remove variables with p-values > 0.2
binom.ra <- glm(class ~ ra, family = binomial, data = SDSS)
summary(binom.ra)
binom.dec <- glm(class ~ dec, family = binomial, data = SDSS)
summary(binom.dec)
binom.u <- glm(class ~ u, family = binomial, data = SDSS)
summary(binom.u)
binom.g <- glm(class ~ g, family = binomial, data = SDSS)
summary(binom.g)
binom.r <- glm(class ~ r, family = binomial, data = SDSS)
summary(binom.r)
binom.i <- glm(class ~ i, family = binomial, data = SDSS)
summary(binom.i)
binom.z <- glm(class ~ z, family = binomial, data = SDSS)
summary(binom.z)
binom.run <- glm(class ~ run, family = binomial, data = SDSS)
summary(binom.run)
binom.field <- glm(class ~ field, family = binomial, data = SDSS)
summary(binom.field)
binom.plate <- glm(class ~ plate, family = binomial, data = SDSS)
summary(binom.plate)
binom.mjd <- glm(class ~ mjd, family = binomial, data = SDSS)
summary(binom.mjd)
binom.fiberid <- glm(class ~ fiberid, family = binomial, data = SDSS)
summary(binom.fiberid)
# not significant: i, field

# initial main effects model
binom.me <- glm(class ~ ra + dec + u + g + r + z + run + camcol + redshift + plate + mjd 
                + fiberid, family = binomial, data = SDSS)
summary(binom.me)

# conduct backward elimination
step(binom.me, direction = "backward")
binom.step <- glm(class ~ dec + r + camcol + redshift, family = binomial, 
                   data = SDSS)
summary(binom.step)

# add additional significant variables back into the model
binom.ra2 <- glm(class ~ ra + dec + r + camcol + redshift, family = binomial, 
                  data = SDSS)
summary(binom.ra2)
binom.u2 <- glm(class ~ u + dec + r + camcol + redshift, family = binomial, 
                 data = SDSS)
summary(binom.u2)
binom.g2 <- glm(class ~ g + dec + r + camcol + redshift, family = binomial, 
                 data = SDSS)
summary(binom.g2)
binom.r2 <- glm(class ~ r + dec + z + camcol + redshift + plate, family = binomial, 
                 data = SDSS)
summary(binom.r2)
binom.i2 <- glm(class ~ i + dec + r + camcol + redshift, family = binomial, 
                 data = SDSS)
summary(binom.i2)
binom.run2 <- glm(class ~ run + dec + r + camcol + redshift, family = binomial, 
                 data = SDSS)
summary(binom.run2)
binom.field2 <- glm(class ~ field + dec + r + camcol + redshift, family = binomial, 
                 data = SDSS)
summary(binom.field2)
binom.mjd2 <- glm(class ~ mjd + dec + r + camcol + redshift, family = binomial, 
                 data = SDSS)
summary(binom.mjd2)
binom.fiberid2 <- glm(class ~ fiberid + dec + r + camcol + redshift, family = binomial, 
                 data = SDSS)
summary(binom.fiberid2)
# none significant

# check for plausible interactions
binom.dr <- glm(class ~ dec + r + camcol + redshift + dec*r, family = binomial, 
                data = SDSS)
summary(binom.dz)
binom.dc <- glm(class ~ dec + r + camcol + redshift + dec*camcol, 
                family = binomial, data = SDSS)
summary(binom.dc)
# significant
binom.dr <- glm(class ~ dec + r + camcol + redshift + dec*redshift, 
                family = binomial, data = SDSS)
summary(binom.dr)
binom.rc <- glm(class ~ dec + r + camcol + redshift + r*camcol, 
                family = binomial, data = SDSS)
summary(binom.rc)
binom.rr <- glm(class ~ dec + r + camcol + redshift + r*redshift, 
                family = binomial, data = SDSS)
summary(binom.rr)
binom.cr <- glm(class ~ dec + r + camcol + redshift + camcol*redshift, 
                family = binomial, data = SDSS)
summary(binom.cr)
binom.cp <- glm(class ~ dec + r + camcol + redshift + camcol*plate, 
                family = binomial, data = SDSS)
summary(binom.cp)

# final model
binom.final <- glm(class ~ dec + r + camcol + redshift + dec*redshift, 
                   family = binomial, data = SDSS)
summary(binom.final)

# check the significance of the final model against a null model
binom.fit0 <- glm(class ~ 1, family = binomial, data = SDSS)
anova(binom.final, binom.fit0, test = "LRT")


# create a confusion matrix and check the accuracy of both models
pred <- predict(fit.final, test, type = "response")
bcl.pred <- data.frame(pred)
bcl.classes <- max.col(bcl.pred) - 1

test[test=="STAR"] <- 0
test[test=="GALAXY"] <- 1
test[test=="QSO"] <- 2

library(cvms)
conf_mat1 <- confusion_matrix(targets = test$class, predictions = bcl.classes)
plot_confusion_matrix(conf_mat1$`Confusion Matrix`[[1]])

library(caret)
confusionMatrix(as.factor(test$class), as.factor(bcl.classes))
# baseline category logit model accuracy: 98.37%

test[test==2] <- 1
binom.pred <- predict(binom.final, test, type = "response")
binom.classes <- ifelse(binom.pred > 0.5, 1, 0)

conf_mat2 <- confusion_matrix(targets = test$class, predictions = binom.classes)
plot_confusion_matrix(conf_mat2$`Confusion Matrix`[[1]])

confusionMatrix(as.factor(test$class), as.factor(binom.classes))
# binomial logistic regression accuracy: 99.77%
