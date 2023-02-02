# plotting 

# load packages
library(sjPlot)
library(sjmisc)
library(sjlabelled)

# load sample data set.
data(efc)
head(efc)
head(sleepstudy)

library(lme4)
fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fit
# prepare group variable
efc$grp = as.factor(efc$e15relat)
levels(x = efc$grp) <- get_labels(efc$e15relat)


# data frame for fitted model
mydf <- data.frame(neg_c_7 = efc$neg_c_7,
                   sex = to_factor(efc$c161sex),
                   c12hour = efc$c12hour,
                   barthel = efc$barthtot,
                   grp = efc$grp)

head(mydf)
unique(mydf$grp)
# fit 2nd model
fit2 <- lmer(neg_c_7 ~ sex + c12hour + barthel + (1| grp), data = mydf)
summary(fit2)
plot_model(fit)

plot_model(fit2, type = "slope", vars = "c12hour")
# ?plot_model

ests <- sjPlot::plot_model(fit2,
                           show.values=TRUE, show.p=TRUE)

ests

plot_model(fit2,type="pred",
                  terms=c("c12hour","grp"))

fm1 <- lmer("Reaction ~ Days + (Days | Subject)", sleepstudy)

sjPlot::plot_model(fm1, type="pred", terms=c("Days","Subject"), pred.type="re", ci.lvl=NA)
