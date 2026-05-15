library(lme4)
library(factoextra)
library(FactoMineR)
library(tidyr)
library(ade4)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(lubridate)
library(MASS)
library(RColorBrewer)
library(wesanderson)
library(MuMIn)
library(ggpubr)
library(sjPlot)
library(nlme)
library(lmerTest)
library(tibble)
library(psych)
library(gridExtra)
library(DHARMa)
library(glmmTMB)
library(AICcmodavg)
library(Factoshiny)
library(FactoInvestigate)
library(installr)
library(aod)
library(maps)

## Données from EcoTaxa
# c'est ce qu'on a obtenu après zooscan + ecotaxa
# toujours repartir d'ici si on veut changer quelque chose
# Ne jamais le faire sur ecotaxa

ecotaxa_brute <- read.csv(file = "Ecotaxa.csv", header = T, sep = ";")
View(ecotaxa_brute)

## Pour ce genre d'étude, seules quelques colonnes sont utiles à prendre en compte parmi les 160 présentes
## On veut garder les colonnes Id, Lake, Communauté, Catégorie, Taille major et Taille féret

ecotaxa <- dplyr::select(ecotaxa_brute,ID,Lake,object_annotation_category,Community,object_major)

## Trie des lignes, selection des groupes de zooplancton seulement, on retire les groupes liés aux doutes 
## ou aux déchets par exemple

ecotaxa <- rbind(
  (Calanoides <- dplyr::filter(ecotaxa,object_annotation_category == "Calanoida")),
  (Cyclopoides <- dplyr::filter(ecotaxa,object_annotation_category == "Cyclopoida")),
  (Daphnies <- dplyr::filter(ecotaxa,object_annotation_category == "Daphnia<Daphniidae")),
  (Bosminidae <- dplyr::filter(ecotaxa,object_annotation_category == "Bosminidae")),
  (Holopediums <- dplyr::filter(ecotaxa,object_annotation_category =="Holopediidae")),
  (Diaphanosoma <- dplyr::filter(ecotaxa,object_annotation_category =="Diaphanosoma")),
  (Chaoborus <- dplyr::filter(ecotaxa,object_annotation_category =="Chaoboridae")),
  (Acari <- dplyr::filter(ecotaxa,object_annotation_category =="Acari")),
  (Sida <- dplyr::filter(ecotaxa,object_annotation_category =="Sida")),
  (Polyphemus <- dplyr::filter(ecotaxa,object_annotation_category =="Polyphemus")),
  (Asplanchna <- dplyr::filter(ecotaxa,object_annotation_category =="Asplanchna")),
  (Ceriodpahnia  <- dplyr::filter(ecotaxa,object_annotation_category =="Ceriodaphnia")),
  (Conochilus <- dplyr::filter(ecotaxa,object_annotation_category =="Conochilus")),
  (Leptodora <- dplyr::filter(ecotaxa,object_annotation_category =="Leptodora kindtii")),
  (Chironomides <- dplyr::filter(ecotaxa,object_annotation_category =="Chironomidae"))
)

## J'ai à la fois trié les lignes du tableau mais aussi créer des dataframes pour chaque groupe de zooplankton

## Maintenant on regroupe les communautés en une seule
## celles avec SAFO only et celles avec au moins une sp compétitrice SEAT CACO

CommunityF <- ifelse(ecotaxa$Community=="SAFO","SAFO","SEAT-CACO")

ecotaxa <- ecotaxa %>%
  dplyr::select(-Community) %>%
  dplyr::mutate(Community = CommunityF)

## On veut que nos valeurs de major soit exprimées en mm et non en pixels
## On sait, via la notice du zooscan, que 1 pixel = 0,0106 mm

ecotaxa <- ecotaxa %>%
  dplyr::mutate(Major_mm = object_major * 0.0106) %>%
  dplyr::select(-object_major)

## Valeur seuil des cyclopoides prédateurs = 800 µm

ecotaxaSeuil <- ecotaxa %>%
  dplyr::mutate(Taille0.8 = case_when(
    object_annotation_category %in% c("Cyclopoida", "Calanoida") & Major_mm <= 0.8 ~ "Petit",
    object_annotation_category %in% c("Cyclopoida", "Calanoida") & Major_mm > 0.8  ~ "Grand",
    TRUE ~ NA_character_
  ))

Cyclo0.8 <- ecotaxaSeuil %>% filter(object_annotation_category == "Cyclopoida")
Cala0.8 <- ecotaxaSeuil %>% filter(object_annotation_category == "Calanoida")

summary(Cyclo0.8$Lake[Cyclo0.8$Taille0.8=="Petit"])
summary(Cyclo0.8$Lake[Cyclo0.8$Taille0.8=="Grand"])

# 713 cyclopoides sont dans la catégories "Grands" quand Vseuil = 800 µm

## Valeur seuil des cyclopoides prédateurs = 1 mm

ecotaxaSeuil <- ecotaxaSeuil %>%
  dplyr::mutate(Taille1 = case_when(
    object_annotation_category %in% c("Cyclopoida", "Calanoida") & Major_mm <= 1 ~ "Petit",
    object_annotation_category %in% c("Cyclopoida", "Calanoida") & Major_mm > 1  ~ "Grand",
    TRUE ~ NA_character_
  ))

Cyclo1<- ecotaxaSeuil %>% filter(object_annotation_category == "Cyclopoida")
Cala1 <- ecotaxaSeuil %>% filter(object_annotation_category == "Calanoida")

summary(Cyclo1$Lake[Cyclo1$Taille1=="Petit"])
summary(Cyclo1$Lake[Cyclo1$Taille1=="Grand"])

# seulment 305 cycloipoides sont "Grands" quand Vseuil = 1 mm 

table(Cyclo0.8$Lake, Cyclo0.8$Taille0.8)
table(Cyclo1$Lake, Cyclo1$Taille1)

write.table(ecotaxaSeuil, "AbZooSeuil.CSV", row.names = FALSE, sep = ";", dec = ",")

# Nous avons nos tableaux zooplancton avec nos variables associées

## Il faut maintenant réunir nos deux tableaux (zoo et parasites)

# Le tableau parasite contient aussi les variables environnementales

# Ici nous allons faire la moyenne des deux réplicats zooplanctons (cob1,cob2)
# Cela revient à considérer que l'unité d'échantillonnage zoo est le lac, 
# pas le filet — ce qui est biologiquement cohérent puisque les poissons
# sont eux aussi échantillonnés à l'échelle du lac.
# Les deux filets deviennent simplement des pseudo-réplicats techniques 
# qui augmentent la précision de ton estimation zoo par lac, 
# pas des unités indépendantes.

# Lire les fichiers
zoo <- read.csv2("AbZooSeuil.CSV", stringsAsFactors = FALSE)
parasites <- read.csv2("Parasites.csv", stringsAsFactors = FALSE)

# Cutaway n'a pas le même nom dans les deux tables
zoo <- zoo %>%
  mutate(Lake = recode(Lake, "Cuttaway" = "Cutaway"))

# Étape 1 : agréger le zoo par lac
zoo_agregee <- zoo %>%
  group_by(Lake, ID) %>%
  summarise(
    Cyclo_total    = sum(object_annotation_category == "Cyclopoida"),
    Cyclo_grand0.8 = sum(object_annotation_category == "Cyclopoida" & Taille0.8 == "Grand"),
    Cyclo_grand1   = sum(object_annotation_category == "Cyclopoida" & Taille1 == "Grand"),
    Cala_total     = sum(object_annotation_category == "Calanoida"),
    .groups = "drop"
  ) %>%
  group_by(Lake) %>%
  summarise(
    Cyclo_total    = mean(Cyclo_total),
    Cyclo_grand0.8 = mean(Cyclo_grand0.8),
    Cyclo_grand1   = mean(Cyclo_grand1),
    Cala_total     = mean(Cala_total),
    .groups = "drop"
  )

# Étape 2 : filtrer Apophallus seulement
apo <- parasites %>%
  filter(Sp_para == "Apophallus")

# Nous allons virer Baie verte qui a des NA
# Et qui n'a pas beaucoup de données d'Apophallus

apo <- apo %>%
  filter(Lake != "Baie_Verte")

# Étape 3 : jointure
data_model <- left_join(apo, zoo_agregee, by = "Lake")

# Étape 4 : variable prévalence
data_model <- data_model %>%
  mutate(Presence = ifelse(Parasites > 0, 1, 0))

# Vérification
sum(is.na(data_model$Cyclo_grand))
setdiff(unique(apo$Lake), unique(zoo$Lake))

# Les deux tableaux sont bien fusionnés

# Le volume d'eau filtré étant le même pour tous les lacs, on est bon comme ça
# Nous avons une abondance relative, censée être représentatif du lac, du zoo
# On va standardiser cette valeur pour faciliter la comparaison des effets 
# dans les modèles sans extrapolation hasardeuse

data_model <- data_model %>%
  mutate(
    Cyclo_total_logstd      = as.numeric(scale(log(Cyclo_total + 1))),
    Cyclo_grand0.8_logstd   = as.numeric(scale(log(Cyclo_grand0.8 + 1))),
    Cyclo_grand1_logstd     = as.numeric(scale(log(Cyclo_grand1 + 1))),
    Cala_total_logstd       = as.numeric(scale(log(Cala_total + 1)))
  )

mean(data_model$Cyclo_grand1_logstd, na.rm = TRUE)
sd(data_model$Cyclo_grand1_logstd, na.rm = TRUE)

data_model <- data_model %>%
  mutate(LF = as.numeric(gsub(",", ".", LF)),
         LF_std = scale(log(LF)))

data_model <- data_model %>%
  mutate(Superficie = as.numeric(gsub(",", ".", Superficie)))

data_model <- data_model %>%
  mutate(Superficie_logstd = scale(log(Superficie)))

data_model <- data_model %>%
  mutate(LF_std = scale(log(LF)))


# Données log transformées et standardisées

## Construction des modèles 

# On met bien toujours LF std en covariable et Lake en effet aléatoire

sum(data_model$Parasites == 0) / nrow(data_model)

abondance0 <- glmer.nb(Parasites ~ LF_std + (1|Lake), data=data_model)
summary(abondance0)
abondance1 <- glmer.nb(Parasites ~ Type + LF_std + (1|Lake), data=data_model)
summary(abondance1)
abondance2 <- glmer.nb(Parasites ~ Superficie_logstd + LF_std + (1|Lake), data=data_model)
summary(abondance2)
abondance3 <- glmer.nb(Parasites ~ Cyclo_grand0.8_logstd + LF_std + (1|Lake), data=data_model)
summary(abondance3)
abondance4 <- glmer.nb(Parasites ~ Cyclo_total_logstd + LF_std + (1|Lake), data=data_model)
summary(abondance4)
abondance5 <- glmer.nb(Parasites ~ Superficie_logstd + Type + LF_std + (1|Lake), data=data_model)
summary(abondance5)
abondance6 <- glmer.nb(Parasites ~ Superficie_logstd + Type + Cyclo_grand0.8_logstd 
+ LF_std +(1|Lake), data=data_model)
summary(abondance6)

# Donc le modèle avec les cyclopoides supérieurs à 800 µm n'améliore pas
# Le modèle le meilleur est celui avec seulement Type/Superficie et biensur LF
# On peut se poser la question si ce n'est pasl a faute du seuil des cyclo
# Moins on va en prendre, plus on va considérer des cyclo grands avec un effet
# probablement plus fort que si on inclut les plus petits

abondance7 <- glmer.nb(Parasites ~ Cyclo_grand1_logstd + LF_std + (1|Lake), data=data_model)
summary(abondance7)
abondance8 <- glmer.nb(Parasites ~ Cyclo_grand1_logstd + Superficie_logstd + LF_std + (1|Lake), data=data_model)
summary(abondance8)
abondance9 <- glmer.nb(Parasites ~ Cyclo_grand1_logstd + Superficie_logstd + Type
+ LF_std + (1|Lake), data=data_model)
summary(abondance9)
abondance10 <- glmer.nb(Parasites ~ Cyclo_grand1_logstd + Type + LF_std + (1|Lake), data=data_model)
summary(abondance10)

cor(data_model$Cyclo_grand1_logstd, data_model$Superficie_logstd, use = "complete.obs")
cor(data_model$Cyclo_grand1_logstd, as.numeric(as.factor(data_model$Type)), use = "complete.obs")

prevalence1 <- glmer(Presence ~ LF_std + (1|Lake), 
                     family = binomial, 
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 2e5)),
                     data = data_model)
summary(prevalence1)

prevalence2 <-glmer(Presence ~ Type + LF_std + (1|Lake), 
                    family = binomial, 
                    data = data_model)
summary(prevalence2)

prevalence3 <- glmer(Presence ~ LF_std + Superficie_logstd + (1|Lake), 
                     family = binomial, 
                     data = data_model)
summary(prevalence3)

prevalence4 <- glmer(Presence ~ Cyclo_grand1_logstd + LF_std + (1|Lake),
                     family = binomial,
                     data = data_model,
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 2e5)))

summary(prevalence4)

prevalence5 <- glmer(Presence ~ Superficie_logstd + Type + LF_std + (1|Lake),
                     family = binomial,
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 2e5)),
                     data = data_model)

summary(prevalence5)

prevalence6 <- glmer(Presence ~ Superficie_logstd + Type + Cyclo_grand1_logstd + LF_std + (1|Lake),
                     family = binomial,
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 2e5)),
                     data = data_model)

summary(prevalence6)

