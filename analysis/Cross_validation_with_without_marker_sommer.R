#########################################
# Cross validation for the trials
##########################################
# split data to training and testing sets
library(tidyverse)
library(sommer)
library(caret)
data("oats")
head(oats)
#taining set - 80% of the data
train.index <- oats$Y %>%
  createDataPartition(p = 0.8, list = FALSE)




train.control <- trainControl(method = "cv", number = 5)
fold5 = caret::createMultiFolds(y = oats$Y, k = 5, times = 100)
is.list(fold5)
for(i in length(fold5)){
index = fold5[[i]] # the index of the sample for training set

train.data <- oats %>%
  filter(row_number() %in% index) # subset the training set
test.data <- oats %>%
  filter(!row_number() %in% index) # subset the testing set

ans1 <- mmer(Y~1,
             random= ~ V + N,
             rcov= ~ units,
             data=train.data)

}

# testing set


# Verify the balance between the training and testing set on the response variable
ggplot()+
  geom_density(data = train.data, aes(x = Y),
               fill = "#00BCD8", alpha = 0.3) +
  geom_density(data = test.data, mapping = aes(x = Y),
               fill = "#F8766D", alpha = 0.3)
head(train.data)
## Fit the model using the training set, and compute the RMSE and MAE
##################################################
# Lets try with sommer
#################################################
library(sommer)

ans1 <- mmer(Y~1,
             random= ~ V + N,
             rcov= ~ units,
             data=train.data)

s = summary(ans1)

### Extract the effects of the factor variables in the model
Veffect = ft$dataWithFitted %>% group_by(V) %>%
  summarise_at(.vars = "V.fitted", .funs =  mean)

Neffect = ft$dataWithFitted %>% group_by(N) %>%
  summarise_at(.vars = "N.fitted", .funs =  mean)
Gmeabn = ans1$Beta["Estimate"]

#############################
###Predict the test data
#############################


test.data %>% add_column(Neffect = NA, Veffect = NA, Gmean = NA,
                         fitted.Y = NA) # Add the columns to add the estimated effects of the fixed variable

##############################################
#  Put the estimated effect of each factor variable in the testdata table
##############################################
for(i in as.vector(Neffect$N)){
  for(j in as.vector(Veffct$V)){
  test.data[test.data$N == i,"Neffect"] = Neffect[Neffect$N == i, "N.fitted"]

  test.data[test.data$V == j,"Veffect"] = Veffect[Veffect$V == j, "V.fitted"]
  test.data[,"Gmean"] = mean(train.data$Y)
  test.data[,"fitted.Y"] = test.data[,"Gmean"] +test.data[,"Neffect"] +
                            test.data[,"Veffect"]
  preictability = cor(test.data$Y, test.data$fitted.Y)
  data.frame(RMSE = RMSE(test.data$fitted.Y, test.data$Y),
             MAE = MAE(test.data$fitted.Y, test.data$Y))
  }
  }
}

?data("oats")
ans1$Beta
ans2 <- lm(Y~ V + N,
             data=train.data)
ans2$coefficients
coef(ans2)
mean(train.data$Y)
predictions <- ans2 %>%
  predict(test.data)
cor(predictions,test.data$Y)
data.frame(RMSE = RMSE(predictions, test.data$Y),
           MAE = MAE(predictions, test.data$Y))
fitted(ans1)
coef(ans1)
