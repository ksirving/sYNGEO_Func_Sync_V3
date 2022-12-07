vignette('mdmr-vignette')

install.packages("MDMR")
library(MDMR)

data(mdmrdata)
head(data)

D[1:10]
D
class(D)
D <- dist(Y.mdmr, method = "manhattan")

mdmr.res <- mdmr(X = X.mdmr, D = D)
summary(mdmr.res)
