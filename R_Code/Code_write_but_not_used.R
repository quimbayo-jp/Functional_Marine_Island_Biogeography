# Other models ------
Model_FRichness <- lmer (qlogis(FRic) ~ Age_Ma_Scale + Age_Ma_Scale_2 +
                           Past_Area_Scale + Past_Area_Scale_2 +
                           #Current_Area_Scale + Current_Area_Scale_2 +
                           Isolation_Scale + Isolation_Scale_2 +
                           (1|Realm), REML=TRUE,
                         data = FDIslands_Factors)
drop1 (Model_FRichness, test = "Chisq")
print (summary (Model_FRichness), correlation=FALSE)
plot (hist(residuals(Model_FRichness)))
shapiro.test(residuals(Model_FRichness))
plot (resid(Model_FRichness))
hist (residuals(Model_FRichness))

Model_FRichness <- glmer (FEve ~ Age_Ma_Scale + #Age_Ma_Scale_2 +
                            Past_Area_Scale + #Past_Area_Scale_2 +
                            Current_Area_Scale + #Current_Area_Scale_2 +
                            Isolation_Scale + #Isolation_Scale_2 +
                            (1|Region), family = binomial,
                          data = FDIslands_Factors, verbose = F)
print(summary(Model_FRichness),correlation=FALSE)
plot (resid(Model_FRichness, type="pearson"))

library(glmmTMB)
m <- glmmTMB (FRic ~ Age_Ma_Scale + Age_Ma_Scale_2 + (1|Realm), 
              data=FDIslands_Factors, family = list (family="beta", link="logit"))
summary (m)
hist (residuals(m))
