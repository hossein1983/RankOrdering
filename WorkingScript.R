library(Hmisc)
library(MASS)
library(pROC)
library(corrplot)
library(gridExtra)
library(tidyverse)

########################################################################################################################################

# read in covariance structure
Sigma = as.matrix(read_csv('Covariance.csv', col_names=F))
N <- 100000

########################################################################################################################################

# set the seed
set.seed(121)

# simulate Xs
df <- mvrnorm(n = N,
              mu = rep(0, 8),
              Sigma =Sigma) %>%
  data.frame() %>% tibble()

corrplot(cor(df %>% select(X1:X5)), type='lower', method = 'number',
         title='Correlation Structure among Underlying Causes')

# df <- df %>%
#   mutate(Z = 1 + 3*X1 + X2 - 2*X3 + 0.5*X4 - 1.5*X5 + 2.5*X6 - 4*X7 - 5*X8 - 3*X1*X3 + 4*X4*X7 + rnorm(N,0,5),
#          Y = ifelse(Z > quantile(Z,.85),1,0))

df <- df %>%
  mutate(Z = -6 + 2*X1 + X2 - 1.5*X3 + 3*X4 - 2*X5 - 4*X3*X5,
         P = 1/(1+exp(-Z)),
         Y = rbinom(N, 1, P))

describe(df$Y)

# data.frame(Z = seq(-50,50,.1)) %>%
#            mutate(P = 1/(1+exp(-Z))) %>%
#   ggplot(aes(Z,P))+
#   geom_line()

########################################################################################################################################

# full model
model_full <- glm(Y ~ X1 + X2 + X3 + X4 + X5 + X3*X5,
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
model_orig <- glm(Y ~ X1 + X2 + X3,
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
       color='',
       y = '')


########################################################################################################################################

# Scenario 1: all probabilities go down 60%
df <- df %>%
  mutate(P1 = P * 0.4,
         Y1 = rbinom(N, 1, P1)) 

df_lift %>%
  select(-Full) %>%
  left_join(df %>%
              group_by(Decile = cut(pred_orig, quantile(pred_orig,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
              summarize(SCN1 = mean(Y1)),
            by = 'Decile') %>%
  pivot_longer(!Decile) %>%
  ggplot(aes(x=Decile,y=value,color=name,group=name))+
  geom_line()+
  labs(title = 'Lift Charts',
       color='',
       y = '')

# ROC curve
df %>%
  with(plot(roc(Y1,pred_orig),print.auc=T))

# loop over scalars
sapply(seq(0.05,5,.05), function(i)
  df %>%
    mutate(P = P * i,
           Y = rbinom(N, 1, P)) %>%
    with(auc(roc(Y,pred_orig),print.auc=T))) %>%
  tibble(AUC = ., Prct_Change = seq(0.05,5,.05)) %>%
  ggplot(aes(Prct_Change,AUC))+
  geom_line()

########################################################################################################################################

# Scenario 2: all probabilities are reduced by a flat bump
df <- df %>%
  mutate(P2 = P - 0.2,
         P2 = ifelse(P2 < 0, 0, P2),
         Y2 = rbinom(N, 1, P2))


df_lift %>%
  select(-Full) %>%
  left_join(df %>%
              group_by(Decile = cut(pred_orig, quantile(pred_orig,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
              summarize(SCN1 = mean(Y2)),
            by = 'Decile') %>%
  pivot_longer(!Decile) %>%
  ggplot(aes(x=Decile,y=value,color=name,group=name))+
  geom_line()+
  labs(title = 'Lift Charts',
       color='',
       y = '')

# ROC curve
df %>%
  with(plot(roc(Y2,pred_orig),print.auc=T))

########################################################################################################################################

# Scenario 3: Effect is present only in certain model deciles
df <- df %>%
   mutate(Decile3 = as.numeric(cut(P, quantile(P,seq(0,1,.1)),include.lowest=T,labels=1:10)),
          P3 = ifelse(Decile3 > 5, (1 / Decile3)*P, P),
          P3 = ifelse(P3 < 0, 0, P3),
          P3 = ifelse(P3 > 1, 1, P3),
          Y3 = rbinom(N, 1, P3))


df_lift %>%
  select(-Full) %>%
  left_join(df %>%
              group_by(Decile = cut(pred_orig, quantile(pred_orig,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
              summarize(SCN1 = mean(Y3)),
            by = 'Decile') %>%
  pivot_longer(!Decile) %>%
  ggplot(aes(x=Decile,y=value,color=name,group=name))+
  geom_line()+
  labs(title = 'Lift Charts',
       color='',
       y = '')

# ROC curve
df %>%
  with(plot(roc(Y3,pred_orig),print.auc=T))


########################################################################################################################################

# Scenario 4: change of relationship with unseen variables
df <- df %>%
  mutate(Z4 = -6 + 2*X1 + X2 - 1.5*X3 - 3*X4 + 2*X5 + 4*X3*X5,
         # Z = -6 + 2*X1 + X2 - 1.5*X3 + 3*X4 - 2*X5 - 4*X3*X5,
         P4 = 1/(1+exp(-Z4)),
         Y4 = rbinom(N, 1, P4))

df_lift %>%
  select(-Full) %>%
  left_join(df %>%
              group_by(Decile = cut(pred_orig, quantile(pred_orig,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
              summarize(SCN1 = mean(Y4)),
            by = 'Decile') %>%
  pivot_longer(!Decile) %>%
  ggplot(aes(x=Decile,y=value,color=name,group=name))+
  geom_line()+
  labs(title = 'Lift Charts',
       color='',
       y = '')

# ROC curve
df %>%
  with(plot(roc(Y4,pred_orig),print.auc=T))

########################################################################################################################################

# Scenario 5: change in importance of unseen variables
df <- df %>%
  mutate(Z5 = -6 + 2*X1 + X2 - 1.5*X3 + 10*X4 - 1*X5 - 9*X3*X5,
         # Z = -6 + 2*X1 + X2 - 1.5*X3 + 3*X4 - 2*X5 - 4*X3*X5,
         P5 = 1/(1+exp(-Z5)),
         Y5 = rbinom(N, 1, P5))

df_lift %>%
  select(-Full) %>%
  left_join(df %>%
              group_by(Decile = cut(pred_orig, quantile(pred_orig,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
              summarize(SCN1 = mean(Y5)),
            by = 'Decile') %>%
  pivot_longer(!Decile) %>%
  ggplot(aes(x=Decile,y=value,color=name,group=name))+
  geom_line()+
  labs(title = 'Lift Charts',
       color='',
       y = '')

# ROC curve
df %>%
  with(plot(roc(Y5,pred_orig),print.auc=T))

########################################################################################################################################

# Scenario 6: change in relationship with seen variables
df <- df %>%
  mutate(Z6 = -6 - 2*X1 - X2 + 1.5*X3 + 3*X4 - 2*X5 - 4*X3*X5,
         # Z = -6 + 2*X1 + X2 - 1.5*X3 + 3*X4 - 2*X5 - 4*X3*X5,
         P6 = 1/(1+exp(-Z6)),
         Y6 = rbinom(N, 1, P6))

df_lift %>%
  select(-Full) %>%
  left_join(df %>%
              group_by(Decile = cut(pred_orig, quantile(pred_orig,seq(0,1,.1)),include.lowest=T,labels=1:10)) %>%
              summarize(SCN1 = mean(Y5)),
            by = 'Decile') %>%
  pivot_longer(!Decile) %>%
  ggplot(aes(x=Decile,y=value,color=name,group=name))+
  geom_line()+
  labs(title = 'Lift Charts',
       color='',
       y = '')

# ROC curve
df %>%
  with(plot(roc(Y5,pred_orig),print.auc=T))

# Update model:
model_updated <- glm(Y6 ~ X1 + X2 + X3,
                  data = df,
                  family = 'binomial')

df <- df %>%
  mutate(pred_updated = predict(model_updated,type='response'))

summary(model_updated)

# ROC curve
plot(roc(df$Y6,df$pred_updated), print.auc=T)
