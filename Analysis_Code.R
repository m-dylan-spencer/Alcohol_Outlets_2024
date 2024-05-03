# Packages that are necessary

library(sf)
library(tidyverse)
library(corrr)
library(tmap)
library(spdep)
library(car)
library(spatialreg)
library(flextable)
library(pastecs)

# This reads in the data

data <- st_read("")
colnames(data)

### Descriptive stats ###

data3 <- data %>%
  select("SUM_IDW_Sh", "IDW_ONPrem", "SUM_IDW_gr", 
         "IDW_Drug_s", "SUM_Liquor","Train_Stop","Major_Rds","IQV2",
         "Res_Mobili","Disadv_Fac", "IDW_Murder", "SUM_IDW_As",
         "SUM_IDW_Ro", "TOTAL_VIOL")

data3 <- data3 %>% 
  st_drop_geometry()

stat.desc(data3)

library("haven")

data <- read_sav("")

model1 <- data %>%
  select("FemaleHead", "Poverty", "Pub_Assist", "Low_Educ", "Unemploy")

#model1 <- model1 %>% 
#  st_drop_geometry()

library(Hmisc)

# Compute correlation matrix with p-values
cor_matrix <- rcorr(as.matrix(model1), type = "spearman")

options(digits = 3, scipen = 999)
print(cor_matrix$r)
print(cor_matrix$P)

### Selecting variables for model 1 ###

model1 <- data %>%
  select("SUM_IDW_Sh", "IDW_ONPrem", "SUM_IDW_gr", 
         "IDW_Drug_s", "SUM_Liquor","Train_Stop","Major_Rds","IQV2",
         "Res_Mobili","Disadv_Fac","LAG_SHOOT")

# Calculating correlations

cor.table <- model1 %>%
  st_drop_geometry() %>%
  correlate()

cor.table

# Histogram of dependent variable

model1 %>%
  ggplot() +
  geom_histogram(aes(x=log(SUM_IDW_Sh))) +
  xlab("Prevalance of Shootings")

# Adding a constant to the response variable

model1$log_SUM_IDW_Sh <- log(model1$SUM_IDW_Sh + 1e-10)

# Multiple regression model

fit.ols.multiple.a <- lm(log_SUM_IDW_Sh ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                         SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                         Res_Mobili + Disadv_Fac, data = model1)

summary(fit.ols.multiple.a)

fit.ols.multiple.b <- lm(SUM_IDW_Sh ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                           SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                           Res_Mobili + Disadv_Fac, data = model1)

summary(fit.ols.multiple.b)

ggplot() + 
  geom_histogram(mapping = aes(x=resid(fit.ols.multiple.a))) +
  xlab("OLS residuals")

ggplot() + 
  geom_histogram(mapping = aes(x=resid(fit.ols.multiple.b))) +
  xlab("OLS residuals")

qqPlot(fit.ols.multiple.a)
qqPlot(fit.ols.multiple.b)

plot(resid(fit.ols.multiple.a))
plot(resid(fit.ols.multiple.b))

tm_shape(model1, unit = "mi") +
  tm_polygons(col = "log_SUM_IDW_Sh", style = "fisher",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Shooting Prevalence, Bronx CBG ",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE, 
            attr.outside = TRUE)

tm_shape(model1, unit = "mi") +
  tm_polygons(col = "SUM_IDW_Sh", style = "fisher",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Shooting Prevalence, Bronx CBG ",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE, 
            attr.outside = TRUE)

model1 <- model1 %>%
  mutate(olsresid.a = resid(fit.ols.multiple.a))

model1 <- model1 %>%
  mutate(olsresid.b = resid(fit.ols.multiple.b))

tm_shape(model1, unit = "mi") +
  tm_polygons(col = "olsresid.a", style = "quantile",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Residuals from linear regression in Bronx CBG",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE,
            attr.outside = TRUE)

tm_shape(model1, unit = "mi") +
  tm_polygons(col = "olsresid.b", style = "quantile",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Residuals from linear regression in Bronx CBG",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE,
            attr.outside = TRUE)

### ESDA

b<-poly2nb(model1, queen=T)
summary(b)
w<-nb2listw(b, style="W", zero.policy = TRUE)

moran.plot(as.numeric(scale(model1$SUM_IDW_Sh)), listw=w, 
           xlab="Standardized Homicide Prevalence", 
           ylab="Neighbors Standardized Homicide Prevalence",
           main=c("Moran Scatterplot for Homicide Prevalence", "in Bronx"))

## Global tests of autocorrelation

moran.test(model1$SUM_IDW_Sh, w)
moran.mc(model1$SUM_IDW_Sh, w, nsim=999)
lm.morantest(fit.ols.multiple.b, w)
geary.test(model1$SUM_IDW_Sh, w)
geary.mc(model1$SUM_IDW_Sh, w, nsim = 999)

## Local autocorrelation

# Local Moran's I

set.seed(1918)
locali<-localmoran_perm(model1$SUM_IDW_Sh, w,
                        nsim = 999) %>%
  as_tibble() %>%
  set_names(c("local_i", "exp_i", "var_i", "z_i", "p_i",
              "p_i_sim", "pi_sim_folded", "skewness", "kurtosis"))

model1 <- model1 %>%
  bind_cols(locali)

model1 <- model1 %>%
  mutate(shootratez =  as.numeric(scale(SUM_IDW_Sh)),
         shootlag = lag.listw(w, SUM_IDW_Sh),
         lisa_cluster = case_when(
           p_i >= 0.05 ~ "Not significant",
           shootratez > 0 & local_i > 0 ~ "High-high",
           shootratez > 0 & local_i < 0 ~ "High-low",
           shootratez < 0 & local_i > 0 ~ "Low-low",
           shootratez < 0 & local_i < 0 ~ "Low-high"
         ))

color_values <- c(`High-high` = "red", 
                  `High-low` = "pink", 
                  `Low-low` = "blue", 
                  `Low-high` = "lightblue", 
                  `Not significant` = "white")

ggplot(model1, aes(x = shootratez, 
                   y = shootlag,
                   fill = lisa_cluster)) + 
  geom_point(color = "black", shape = 21, size = 2) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = color_values) + 
  labs(x = "Shooting rate (z-score)",
       y = "Spatial lag of shooting rate (z-score)",
       fill = "Cluster type")


tm_shape(model1, unit = "mi") +
  tm_polygons(col = "lisa_cluster", title = "Local Moran's I", palette = color_values) +
  tm_compass(type = "4star", position = c("right", "bottom")) + 
  tm_scale_bar(breaks = c(0, 1.5, 3), text.size = 1) +
  tm_layout(frame = F, main.title = "Bronx shooting clusters",
            legend.outside = T) 

# G*

localg.self <- include.self(b)
w.self <- nb2listw(localg.self, style = "W")
localgstar <- localG(model1$SUM_IDW_Sh, w.self)
model1 <- model1 %>%
  mutate(localgstar = as.numeric(localgstar))

model1 <- model1 %>%
  mutate(hotspotsgs = case_when(
    localgstar <= -2.58 ~ "Cold spot 99%",
    localgstar > -2.58 & localgstar <=-1.96 ~ "Cold spot 95%",
    localgstar > -1.96 & localgstar <= -1.65 ~ "Cold spot 90%",
    localgstar >= 1.65 & localgstar < 1.96 ~ "Hot spot 90%",
    localgstar >= 1.96 & localgstar <= 2.58 ~ "Hot spot 95%",
    localgstar >= 2.58 ~ "Hot spot 99%",
    TRUE ~ "Not Significant"),
    hotspotsgs = factor(hotspotsgs,
                        levels = c("Cold spot 99%", "Cold spot 95%",
                                   "Cold spot 90%", "Not Significant",
                                   "Hot spot 90%", "Hot spot 95%",
                                   "Hot spot 99%")))

tm_shape(model1, unit = "mi") +
  tm_polygons(col = "hotspotsgs", title = "Gi* value", palette = c("blue","white", "red")) +
  tm_compass(type = "4star", position = c("right", "bottom")) + 
  tm_layout(frame = F, main.title = "Bronx shooting clusters",
            legend.outside = T)

### Spatial lag regression

fit.lag<-lagsarlm(SUM_IDW_Sh ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                    SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                    Res_Mobili + Disadv_Fac, data = model1, 
                  listw = w)
summary(fit.lag)
fit.lag.effects <- impacts(fit.lag, listw = w, R = 999)
options(digits = 3, scipen = 999)
fit.lag.effects
summary(fit.lag.effects, zstats = TRUE, short = T)

### Spatial error model

fit.err <- errorsarlm(SUM_IDW_Sh ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                        SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                        Res_Mobili + Disadv_Fac, data = model1, 
                      listw = w)
summary(fit.err)

### Picking a model

AIC(fit.ols.multiple.b)
AIC(fit.lag)
AIC(fit.err)

# Save AIC values

AICs<-c(AIC(fit.ols.multiple.b),AIC(fit.lag), AIC(fit.err))
labels<-c("OLS","SLM", "SEM")

flextable(data.frame(Models=labels, AIC=round(AICs, 2)))

# Lagrange multiplier test

lm.LMtests(fit.ols.multiple.b, listw = w, test = "all",  zero.policy=TRUE)

### Presenting results

summary(fit.ols.multiple.b)

fit.ols.multiple.b %>%
  as_flextable()

### Model 2 ###

model2 <- data %>%
  select("IDW_Murder", "IDW_ONPrem", "SUM_IDW_gr", 
         "IDW_Drug_s", "SUM_Liquor","Train_Stop","Major_Rds","IQV2",
         "Res_Mobili","Disadv_Fac","LAG_HOMIC")

# Calculating correlations

cor.table <- model2 %>%
  st_drop_geometry() %>%
  correlate()

cor.table

# Histogram of dependent variable

model2 %>%
  ggplot() +
  geom_histogram(aes(x=log(IDW_Murder))) +
  xlab("Prevalance of Homicides")

# Multiple regression model

fit.ols.multiple <- lm(IDW_Murder ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                         SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                         Res_Mobili + Disadv_Fac, data = model2)

summary(fit.ols.multiple)

ggplot() + 
  geom_histogram(mapping = aes(x=resid(fit.ols.multiple))) +
  xlab("OLS residuals")

qqPlot(fit.ols.multiple)

plot(resid(fit.ols.multiple))

tm_shape(model2, unit = "mi") +
  tm_polygons(col = "IDW_Murder", style = "fisher",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Homicide Prevalence, Bronx CBG ",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE, 
            attr.outside = TRUE)

model2 <- model2 %>%
  mutate(olsresid = resid(fit.ols.multiple))

tm_shape(model2, unit = "mi") +
  tm_polygons(col = "olsresid", style = "quantile",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Residuals from linear regression in Bronx CBG",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE,
            attr.outside = TRUE)

### ESDA

b<-poly2nb(model2, queen=T)
summary(b)
w<-nb2listw(b, style="W", zero.policy = TRUE)

moran.plot(as.numeric(scale(model2$IDW_Murder)), listw=w, 
           xlab="Standardized Homicide Prevalence", 
           ylab="Neighbors Standardized Homicide Prevalence",
           main=c("Moran Scatterplot for Homicide Prevalence", "in Bronx"))

## Global tests of autocorrelation

moran.test(model2$IDW_Murder, w)
moran.mc(model2$IDW_Murder, w, nsim=999)
lm.morantest(fit.ols.multiple, w)
geary.test(model2$IDW_Murder, w)
geary.mc(model2$IDW_Murder, w, nsim = 999)

## Local autocorrelation

# Local Moran's I

set.seed(1918)
locali<-localmoran_perm(model2$IDW_Murder, w,
                        nsim = 999) %>%
  as_tibble() %>%
  set_names(c("local_i", "exp_i", "var_i", "z_i", "p_i",
              "p_i_sim", "pi_sim_folded", "skewness", "kurtosis"))

model2 <- model2 %>%
  bind_cols(locali)

model2 <- model2 %>%
  mutate(murderratez =  as.numeric(scale(IDW_Murder)),
         murderlag = lag.listw(w, IDW_Murder),
         lisa_cluster = case_when(
           p_i >= 0.05 ~ "Not significant",
           murderratez > 0 & local_i > 0 ~ "High-high",
           murderratez > 0 & local_i < 0 ~ "High-low",
           murderratez < 0 & local_i > 0 ~ "Low-low",
           murderratez < 0 & local_i < 0 ~ "Low-high"
         ))

color_values <- c(`High-high` = "red", 
                  `High-low` = "pink", 
                  `Low-low` = "blue", 
                  `Low-high` = "lightblue", 
                  `Not significant` = "white")

ggplot(model2, aes(x = murderratez, 
                       y = murderlag,
                       fill = lisa_cluster)) + 
  geom_point(color = "black", shape = 21, size = 2) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = color_values) + 
  labs(x = "Homicide rate (z-score)",
       y = "Spatial lag of homicide rate (z-score)",
       fill = "Cluster type")


tm_shape(model2, unit = "mi") +
  tm_polygons(col = "lisa_cluster", title = "Local Moran's I", palette = color_values) +
  tm_compass(type = "4star", position = c("left", "bottom")) + 
  tm_scale_bar(breaks = c(0, 1.5, 3), text.size = 1) +
  tm_layout(frame = F, main.title = "Bronx homicide clusters",
            legend.outside = T) 

# G*

localg.self <- include.self(b)
w.self <- nb2listw(localg.self, style = "W")
localgstar <- localG(model2$IDW_Murder, w.self)
model2 <- model2 %>%
  mutate(localgstar = as.numeric(localgstar))

model2 <- model2 %>%
  mutate(hotspotsgs = case_when(
    localgstar <= -2.58 ~ "Cold spot 99%",
    localgstar > -2.58 & localgstar <=-1.96 ~ "Cold spot 95%",
    localgstar > -1.96 & localgstar <= -1.65 ~ "Cold spot 90%",
    localgstar >= 1.65 & localgstar < 1.96 ~ "Hot spot 90%",
    localgstar >= 1.96 & localgstar <= 2.58 ~ "Hot spot 95%",
    localgstar >= 2.58 ~ "Hot spot 99%",
    TRUE ~ "Not Significant"),
    hotspotsgs = factor(hotspotsgs,
                        levels = c("Cold spot 99%", "Cold spot 95%",
                                   "Cold spot 90%", "Not Significant",
                                   "Hot spot 90%", "Hot spot 95%",
                                   "Hot spot 99%")))

tm_shape(model2, unit = "mi") +
  tm_polygons(col = "hotspotsgs", title = "Gi* value", palette = c("blue","white", "red")) +
  tm_compass(type = "4star", position = c("right", "bottom")) + 
  tm_layout(frame = F, main.title = "Bronx homicide clusters",
            legend.outside = T) 


### Spatial lag regression

fit.lag<-lagsarlm(IDW_Murder ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                    SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                    Res_Mobili + Disadv_Fac, data = model2, 
                  listw = w)
summary(fit.lag)
options(digits = 3, scipen = 999)
fit.lag.effects <- impacts(fit.lag, listw = w, R = 999)
fit.lag.effects
summary(fit.lag.effects, zstats = TRUE, short = TRUE)

### Spatial error model

fit.err <- errorsarlm(IDW_Murder ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                        SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                        Res_Mobili + Disadv_Fac, data = model2, 
                      listw = w)
summary(fit.err)

### Picking a model

AIC(fit.ols.multiple)
AIC(fit.lag)
AIC(fit.err)

# Save AIC values

AICs<-c(AIC(fit.ols.multiple),AIC(fit.lag), AIC(fit.err))
labels<-c("OLS","SLM", "SEM")

flextable(data.frame(Models=labels, AIC=round(AICs, 2)))

# Lagrange multiplier test

lm.LMtests(fit.ols.multiple, listw = w, test = "all",  zero.policy=TRUE)

### Presenting results

summary(fit.ols.multiple)

fit.ols.multiple %>%
  as_flextable()

### Model 3 ###

model3 <- data %>%
  select("SUM_IDW_As", "IDW_ONPrem", "SUM_IDW_gr", 
         "IDW_Drug_s", "SUM_Liquor","Train_Stop","Major_Rds","IQV2",
         "Res_Mobili","Disadv_Fac","LAG_ASSAUL")

# Calculating correlations

cor.table <- model3 %>%
  st_drop_geometry() %>%
  correlate()

cor.table

# Histogram of dependent variable

model3 %>%
  ggplot() +
  geom_histogram(aes(x=log(SUM_IDW_As))) +
  xlab("Prevalance of Assault")

# Multiple regression model

fit.ols.multiple <- lm(SUM_IDW_As ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                         SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                         Res_Mobili + Disadv_Fac + LAG_ASSAUL, data = model3)

summary(fit.ols.multiple)

ggplot() + 
  geom_histogram(mapping = aes(x=resid(fit.ols.multiple))) +
  xlab("OLS residuals")

qqPlot(fit.ols.multiple)

plot(resid(fit.ols.multiple))

tm_shape(model3, unit = "mi") +
  tm_polygons(col = "SUM_IDW_As", style = "fisher",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Assault Prevalence, Bronx CBG ",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE, 
            attr.outside = TRUE)

model3 <- model3 %>%
  mutate(olsresid = resid(fit.ols.multiple))

tm_shape(model3, unit = "mi") +
  tm_polygons(col = "olsresid", style = "quantile",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Residuals from linear regression in Bronx CBG",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE,
            attr.outside = TRUE)

### ESDA

b<-poly2nb(model3, queen=T)
summary(b)
w<-nb2listw(b, style="W", zero.policy = TRUE)

moran.plot(as.numeric(scale(model3$SUM_IDW_As)), listw=w, 
           xlab="Standardized Assault Prevalence", 
           ylab="Neighbors Standardized Assault Prevalence",
           main=c("Moran Scatterplot for Assault Prevalence", "in Bronx"))

## Global tests of autocorrelation

moran.test(model3$SUM_IDW_As, w)
moran.mc(model3$SUM_IDW_As, w, nsim=999)
lm.morantest(fit.ols.multiple, w)
geary.test(model3$SUM_IDW_As, w)
geary.mc(model3$SUM_IDW_As, w, nsim = 999)

## Local autocorrelation

# Local Moran's I

set.seed(1918)
locali<-localmoran_perm(model3$SUM_IDW_As, w,
                        nsim = 999) %>%
  as_tibble() %>%
  set_names(c("local_i", "exp_i", "var_i", "z_i", "p_i",
              "p_i_sim", "pi_sim_folded", "skewness", "kurtosis"))

model3 <- model3 %>%
  bind_cols(locali)

model3 <- model3 %>%
  mutate(assratez =  as.numeric(scale(SUM_IDW_As)),
         asslag = lag.listw(w, SUM_IDW_As),
         lisa_cluster = case_when(
           p_i >= 0.05 ~ "Not significant",
           assratez > 0 & local_i > 0 ~ "High-high",
           assratez > 0 & local_i < 0 ~ "High-low",
           assratez < 0 & local_i > 0 ~ "Low-low",
           assratez < 0 & local_i < 0 ~ "Low-high"
         ))

color_values <- c(`High-high` = "red", 
                  `High-low` = "pink", 
                  `Low-low` = "blue", 
                  `Low-high` = "lightblue", 
                  `Not significant` = "white")

ggplot(model3, aes(x = assratez, 
                   y = asslag,
                   fill = lisa_cluster)) + 
  geom_point(color = "black", shape = 21, size = 2) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = color_values) + 
  labs(x = "Assault rate (z-score)",
       y = "Spatial lag of assault rate (z-score)",
       fill = "Cluster type")


tm_shape(model3, unit = "mi") +
  tm_polygons(col = "lisa_cluster", title = "Local Moran's I", palette = color_values) +
  tm_compass(type = "4star", position = c("left", "bottom")) + 
  tm_scale_bar(breaks = c(0, 1.5, 3), text.size = 1) +
  tm_layout(frame = F, main.title = "Bronx assault clusters",
            legend.outside = T) 

# G*

localg.self <- include.self(b)
w.self <- nb2listw(localg.self, style = "W")
localgstar <- localG(model3$SUM_IDW_As, w.self)
model3 <- model3 %>%
  mutate(localgstar = as.numeric(localgstar))

model3 <- model3 %>%
  mutate(hotspotsgs = case_when(
    localgstar <= -2.58 ~ "Cold spot 99%",
    localgstar > -2.58 & localgstar <=-1.96 ~ "Cold spot 95%",
    localgstar > -1.96 & localgstar <= -1.65 ~ "Cold spot 90%",
    localgstar >= 1.65 & localgstar < 1.96 ~ "Hot spot 90%",
    localgstar >= 1.96 & localgstar <= 2.58 ~ "Hot spot 95%",
    localgstar >= 2.58 ~ "Hot spot 99%",
    TRUE ~ "Not Significant"),
    hotspotsgs = factor(hotspotsgs,
                        levels = c("Cold spot 99%", "Cold spot 95%",
                                   "Cold spot 90%", "Not Significant",
                                   "Hot spot 90%", "Hot spot 95%",
                                   "Hot spot 99%")))

tm_shape(model3, unit = "mi") +
  tm_polygons(col = "hotspotsgs", title = "Gi* value", palette = c("blue","white", "red")) +
  tm_compass(type = "4star", position = c("right", "bottom")) + 
  tm_layout(frame = F, main.title = "Bronx assault clusters",
            legend.outside = T) 


### Spatial lag regression

fit.lag<-lagsarlm(SUM_IDW_As ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                    SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                    Res_Mobili + Disadv_Fac, data = model3, 
                  listw = w)
summary(fit.lag)
options(digits = 3, scipen = 999)
fit.lag.effects <- impacts(fit.lag, listw = w, R = 999)
fit.lag.effects
summary(fit.lag.effects, zstats = TRUE, short = TRUE)

### Spatial error model

fit.err <- errorsarlm(SUM_IDW_As ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                        SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                        Res_Mobili + Disadv_Fac, data = model3, 
                      listw = w)
summary(fit.err)

### Picking a model

AIC(fit.ols.multiple)
AIC(fit.lag)
AIC(fit.err)

# Save AIC values

AICs<-c(AIC(fit.ols.multiple),AIC(fit.lag), AIC(fit.err))
labels<-c("OLS","SLM", "SEM")

flextable(data.frame(Models=labels, AIC=round(AICs, 2)))

# Lagrange multiplier test

lm.LMtests(fit.ols.multiple, listw = w, test = "all",  zero.policy=TRUE)

### Presenting results

summary(fit.ols.multiple)

fit.ols.multiple %>%
  as_flextable()

### Model 4 ###

model4 <- data %>%
  select("SUM_IDW_Ro", "IDW_ONPrem", "SUM_IDW_gr", 
         "IDW_Drug_s", "SUM_Liquor","Train_Stop","Major_Rds","IQV2",
         "Res_Mobili","Disadv_Fac","LAG_ROBBER")

# Calculating correlations

cor.table <- model4 %>%
  st_drop_geometry() %>%
  correlate()

cor.table

# Histogram of dependent variable

model4 %>%
  ggplot() +
  geom_histogram(aes(x=log(SUM_IDW_Ro))) +
  xlab("Prevalance of Robbery")

# Multiple regression model

fit.ols.multiple <- lm(SUM_IDW_Ro ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                         SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                         Res_Mobili + Disadv_Fac + LAG_ROBBER, data = model4)

summary(fit.ols.multiple)

ggplot() + 
  geom_histogram(mapping = aes(x=resid(fit.ols.multiple))) +
  xlab("OLS residuals")

qqPlot(fit.ols.multiple)

plot(resid(fit.ols.multiple))

tm_shape(model4, unit = "mi") +
  tm_polygons(col = "SUM_IDW_Ro", style = "fisher",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Robbery Prevalence, Bronx CBG ",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE, 
            attr.outside = TRUE)

model4 <- model4 %>%
  mutate(olsresid = resid(fit.ols.multiple))

tm_shape(model4, unit = "mi") +
  tm_polygons(col = "olsresid", style = "quantile",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Residuals from linear regression in Bronx CBG",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE,
            attr.outside = TRUE)

### ESDA

b<-poly2nb(model4, queen=T)
summary(b)
w<-nb2listw(b, style="W", zero.policy = TRUE)

moran.plot(as.numeric(scale(model4$SUM_IDW_Ro)), listw=w, 
           xlab="Standardized Robbery Prevalence", 
           ylab="Neighbors Standardized Robbery Prevalence",
           main=c("Moran Scatterplot for Robbery Prevalence", "in Bronx"))

## Global tests of autocorrelation

moran.test(model4$SUM_IDW_Ro, w)
moran.mc(model4$SUM_IDW_Ro, w, nsim=999)
lm.morantest(fit.ols.multiple, w)
geary.test(model4$SUM_IDW_Ro, w)
geary.mc(model4$SUM_IDW_Ro, w, nsim = 999)

## Local autocorrelation

# Local Moran's I

set.seed(1918)
locali<-localmoran_perm(model4$SUM_IDW_Ro, w,
                        nsim = 999) %>%
  as_tibble() %>%
  set_names(c("local_i", "exp_i", "var_i", "z_i", "p_i",
              "p_i_sim", "pi_sim_folded", "skewness", "kurtosis"))

model4 <- model4 %>%
  bind_cols(locali)

model4 <- model4 %>%
  mutate(robratez =  as.numeric(scale(SUM_IDW_Ro)),
         roblag = lag.listw(w, SUM_IDW_Ro),
         lisa_cluster = case_when(
           p_i >= 0.05 ~ "Not significant",
           robratez > 0 & local_i > 0 ~ "High-high",
           robratez > 0 & local_i < 0 ~ "High-low",
           robratez < 0 & local_i > 0 ~ "Low-low",
           robratez < 0 & local_i < 0 ~ "Low-high"
         ))

color_values <- c(`High-high` = "red", 
                  `High-low` = "pink", 
                  `Low-low` = "blue", 
                  `Low-high` = "lightblue", 
                  `Not significant` = "white")

ggplot(model4, aes(x = robratez, 
                   y = roblag,
                   fill = lisa_cluster)) + 
  geom_point(color = "black", shape = 21, size = 2) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = color_values) + 
  labs(x = "Robbery rate (z-score)",
       y = "Spatial lag of robbery rate (z-score)",
       fill = "Cluster type")

tm_shape(model4, unit = "mi") +
  tm_polygons(col = "lisa_cluster", title = "Local Moran's I", palette = color_values) +
  tm_compass(type = "4star", position = c("left", "bottom")) + 
  tm_scale_bar(breaks = c(0, 1.5, 3), text.size = 1) +
  tm_layout(frame = F, main.title = "Bronx robbery clusters",
            legend.outside = T) 

# G*

localg.self <- include.self(b)
w.self <- nb2listw(localg.self, style = "W")
localgstar <- localG(model4$SUM_IDW_Ro, w.self)
model4 <- model4 %>%
  mutate(localgstar = as.numeric(localgstar))

model4 <- model4 %>%
  mutate(hotspotsgs = case_when(
    localgstar <= -2.58 ~ "Cold spot 99%",
    localgstar > -2.58 & localgstar <=-1.96 ~ "Cold spot 95%",
    localgstar > -1.96 & localgstar <= -1.65 ~ "Cold spot 90%",
    localgstar >= 1.65 & localgstar < 1.96 ~ "Hot spot 90%",
    localgstar >= 1.96 & localgstar <= 2.58 ~ "Hot spot 95%",
    localgstar >= 2.58 ~ "Hot spot 99%",
    TRUE ~ "Not Significant"),
    hotspotsgs = factor(hotspotsgs,
                        levels = c("Cold spot 99%", "Cold spot 95%",
                                   "Cold spot 90%", "Not Significant",
                                   "Hot spot 90%", "Hot spot 95%",
                                   "Hot spot 99%")))

tm_shape(model4, unit = "mi") +
  tm_polygons(col = "hotspotsgs", title = "Gi* value", palette = c("blue","white", "red")) +
  tm_compass(type = "4star", position = c("right", "bottom")) + 
  tm_layout(frame = F, main.title = "Bronx robbery clusters",
            legend.outside = T) 


### Spatial lag regression

fit.lag<-lagsarlm(SUM_IDW_Ro ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                    SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                    Res_Mobili + Disadv_Fac, data = model4, 
                  listw = w)
summary(fit.lag)
options(digits = 2, scipen = 99)
fit.lag.effects <- impacts(fit.lag, listw = w, R = 999)
fit.lag.effects
summary(fit.lag.effects, zstats = TRUE, short = TRUE)

### Spatial error model

fit.err <- errorsarlm(SUM_IDW_Ro ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                        SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                        Res_Mobili + Disadv_Fac, data = model4, 
                      listw = w)
summary(fit.err)

### Picking a model

AIC(fit.ols.multiple)
AIC(fit.lag)
AIC(fit.err)

# Save AIC values

AICs<-c(AIC(fit.ols.multiple),AIC(fit.lag), AIC(fit.err))
labels<-c("OLS","SLM", "SEM")

flextable(data.frame(Models=labels, AIC=round(AICs, 2)))

# Lagrange multiplier test

lm.LMtests(fit.ols.multiple, listw = w, test = "all",  zero.policy=TRUE)

### Presenting results

summary(fit.ols.multiple)

fit.ols.multiple %>%
  as_flextable()

### Model 5 ###

model5 <- data %>%
  select("TOTAL_VIOL", "IDW_ONPrem", "SUM_IDW_gr", 
         "IDW_Drug_s", "SUM_Liquor","Train_Stop","Major_Rds","IQV2",
         "Res_Mobili","Disadv_Fac","LAG_TOTAL")

# Calculating correlations

cor.table <- model5 %>%
  st_drop_geometry() %>%
  correlate()

cor.table

# Histogram of dependent variable

model5 %>%
  ggplot() +
  geom_histogram(aes(x=log(TOTAL_VIOL))) +
  xlab("Prevalance of Violence")

# Multiple regression model

fit.ols.multiple <- lm(TOTAL_VIOL ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                         SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                         Res_Mobili + Disadv_Fac + LAG_TOTAL, data = model5)

summary(fit.ols.multiple)

ggplot() + 
  geom_histogram(mapping = aes(x=resid(fit.ols.multiple))) +
  xlab("OLS residuals")

qqPlot(fit.ols.multiple)

plot(resid(fit.ols.multiple))

tm_shape(model5, unit = "mi") +
  tm_polygons(col = "TOTAL_VIOL", style = "fisher",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Violence Prevalence, Bronx CBG ",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE, 
            attr.outside = TRUE)

model5 <- model5 %>%
  mutate(olsresid = resid(fit.ols.multiple))

tm_shape(model5, unit = "mi") +
  tm_polygons(col = "olsresid", style = "quantile",palette = "Reds", 
              border.alpha = 0, title = "") +
  tm_scale_bar(breaks = c(0, 2, 4), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Residuals from linear regression in Bronx CBG",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE,
            attr.outside = TRUE)

### ESDA

b<-poly2nb(model5, queen=T)
summary(b)
w<-nb2listw(b, style="W", zero.policy = TRUE)

moran.plot(as.numeric(scale(model5$TOTAL_VIOL)), listw=w, 
           xlab="Standardized Violence Prevalence", 
           ylab="Neighbors Standardized Violence Prevalence",
           main=c("Moran Scatterplot for Violence Prevalence", "in Bronx"))

## Global tests of autocorrelation

moran.test(model5$TOTAL_VIOL, w)
moran.mc(model5$TOTAL_VIOL, w, nsim=999)
lm.morantest(fit.ols.multiple, w)
geary.test(model5$TOTAL_VIOL, w)
geary.mc(model5$TOTAL_VIOL, w, nsim = 999)

## Local autocorrelation

# Local Moran's I

set.seed(1918)
locali<-localmoran_perm(model5$TOTAL_VIOL, w,
                        nsim = 999) %>%
  as_tibble() %>%
  set_names(c("local_i", "exp_i", "var_i", "z_i", "p_i",
              "p_i_sim", "pi_sim_folded", "skewness", "kurtosis"))

model5 <- model5 %>%
  bind_cols(locali)

model5 <- model5 %>%
  mutate(vratez =  as.numeric(scale(TOTAL_VIOL)),
         vlag = lag.listw(w, TOTAL_VIOL),
         lisa_cluster = case_when(
           p_i >= 0.05 ~ "Not significant",
           vratez > 0 & local_i > 0 ~ "High-high",
           vratez > 0 & local_i < 0 ~ "High-low",
           vratez < 0 & local_i > 0 ~ "Low-low",
           vratez < 0 & local_i < 0 ~ "Low-high"
         ))

color_values <- c(`High-high` = "red", 
                  `High-low` = "pink", 
                  `Low-low` = "blue", 
                  `Low-high` = "lightblue", 
                  `Not significant` = "white")

ggplot(model5, aes(x = vratez, 
                   y = vlag,
                   fill = lisa_cluster)) + 
  geom_point(color = "black", shape = 21, size = 2) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = color_values) + 
  labs(x = "Violence rate (z-score)",
       y = "Spatial lag of violence rate (z-score)",
       fill = "Cluster type")

tm_shape(model5, unit = "mi") +
  tm_polygons(col = "lisa_cluster", title = "Local Moran's I", palette = color_values) +
  tm_compass(type = "4star", position = c("left", "bottom")) + 
  tm_scale_bar(breaks = c(0, 1.5, 3), text.size = 1) +
  tm_layout(frame = F, main.title = "Bronx violence clusters",
            legend.outside = T) 

# G*

localg.self <- include.self(b)
w.self <- nb2listw(localg.self, style = "W")
localgstar <- localG(model5$TOTAL_VIOL, w.self)
model5 <- model5 %>%
  mutate(localgstar = as.numeric(localgstar))

model5 <- model5 %>%
  mutate(hotspotsgs = case_when(
    localgstar <= -2.58 ~ "Cold spot 99%",
    localgstar > -2.58 & localgstar <=-1.96 ~ "Cold spot 95%",
    localgstar > -1.96 & localgstar <= -1.65 ~ "Cold spot 90%",
    localgstar >= 1.65 & localgstar < 1.96 ~ "Hot spot 90%",
    localgstar >= 1.96 & localgstar <= 2.58 ~ "Hot spot 95%",
    localgstar >= 2.58 ~ "Hot spot 99%",
    TRUE ~ "Not Significant"),
    hotspotsgs = factor(hotspotsgs,
                        levels = c("Cold spot 99%", "Cold spot 95%",
                                   "Cold spot 90%", "Not Significant",
                                   "Hot spot 90%", "Hot spot 95%",
                                   "Hot spot 99%")))

tm_shape(model5, unit = "mi") +
  tm_polygons(col = "hotspotsgs", title = "Gi* value", palette = c("blue","white", "red")) +
  tm_compass(type = "4star", position = c("right", "bottom")) + 
  tm_layout(frame = F, main.title = "Bronx violence clusters",
            legend.outside = T) 


### Spatial lag regression

fit.lag<-lagsarlm(TOTAL_VIOL ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                    SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                    Res_Mobili + Disadv_Fac, data = model5, 
                  listw = w)
summary(fit.lag)
options(digits = 3, scipen = 99)
fit.lag.effects <- impacts(fit.lag, listw = w, R = 999)
fit.lag.effects
summary(fit.lag.effects, zstats = TRUE, short = TRUE)

### Spatial error model

fit.err <- errorsarlm(TOTAL_VIOL ~ IDW_ONPrem + SUM_IDW_gr + IDW_Drug_s +
                        SUM_Liquor + Train_Stop + Major_Rds + IQV2 +
                        Res_Mobili + Disadv_Fac, data = model5, 
                      listw = w)
summary(fit.err)

### Picking a model

AIC(fit.ols.multiple)
AIC(fit.lag)
AIC(fit.err)

# Save AIC values

AICs<-c(AIC(fit.ols.multiple),AIC(fit.lag), AIC(fit.err))
labels<-c("OLS","SLM", "SEM")

flextable(data.frame(Models=labels, AIC=round(AICs, 2)))

# Lagrange multiplier test

lm.LMtests(fit.ols.multiple, listw = w, test = "all",  zero.policy=TRUE)

### Presenting results

summary(fit.ols.multiple)

fit.ols.multiple %>%
  as_flextable()

