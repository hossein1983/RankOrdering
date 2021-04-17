library(Hmisc)
library(MASS)
library(pROC)
library(corrplot)
library(gridExtra)
library(tidyverse)

source('auc_lift.R')

########################################################################################################################################

# read in covariance structure
Sigma = as.matrix(read_csv('Covariance.csv', col_names=F))
N <- 10000

########################################################################################################################################

# set the seed
set.seed(121)

# simulate Xs
df <- mvrnorm(n = N,
              mu = rep(0, 8),
              Sigma =Sigma) %>%
  data.frame() %>% tibble()

corrplot(cor(df %>% select(X1:X8)), type='lower', method = 'number',
         title='Correlation Structure among Underlying Causes')

df <- df %>%
  mutate(Z = 1 - 3*X1 + X2 - 2*X3 + 0.5*X4 - 1.5*X5 + 2.5*X6 - 4**X7 + 5*X8 - 3*X1*X3 + 4*X4*X7 + rnorm(N,0,5),
         Y = ifelse(Z > quantile(Z,.85),1,0))

describe(df$Y)

########################################################################################################################################

# full model
model_full <- glm(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X1*X3 + X4*X7,
                  data = df,
                  family = 'binomial')

df <- df %>%
  mutate(pred_full = predict(model_full,type='response'))

summary(model_full)

# ROC curve
plot(roc(df$Y,df$pred_full), print.auc=T)

# lift data
df_lift <- df %>%
    group_by(Decile = cut(pred_full, quantile(pred_full,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
    summarize(Full = mean(Y))

########################################################################################################################################

# Original model
model_orig <- glm(Y ~ X1 + X2 + X3 + X4 + X5 + X1*X3,
                  data = df,
                  family = 'binomial')

df <- df %>%
  mutate(pred_orig = predict(model_orig,type='response'))

summary(model_orig)

# ROC curve
plot(roc(df$Y,df$pred_orig), print.auc=T)

# add lift data
df_lift <- df_lift %>%
  left_join(df %>%
              group_by(Decile = cut(pred_orig, quantile(pred_orig,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
              summarize(Original = mean(Y)),
            by = 'Decile')

# lift chart
df_lift %>%
  pivot_longer(Full:Original) %>%
  ggplot(aes(x=Decile,y=value,color=name,group=name))+
  geom_line()+
  labs(title = 'Lift Chart for Full and Original Models',
       color='')

