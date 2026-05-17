library(lme4)
library(factoextra)
library(FactoMineR)
library(tidyr)
library(ade4)
library(ggplot2)
library(ggeffects)
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
library(ggmap)
library(maptiles)
library(tidyterra)
library(cowplot)
library(ggrepel)
library(glmmTMB)

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

# Nous allons virer Baie verte qui a des NA pour superficie
# (peut être une idée à la con me direz vous)
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

sum(data_model$Parasites == 0) / nrow(data_model) # donc pas 0 inflated

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

## On va quand même vérfier un peu nos modèles 

sim_res <- simulateResiduals(abondance7)
plot(sim_res)
# pas top

# Modèle zero-inflated
abondance7_zi <- glmmTMB(Parasites ~ Cyclo_grand1_logstd + LF_std + (1|Lake),
                         family = nbinom2,
                         ziformula = ~1,
                         data = data_model)

sim_res_zi <- simulateResiduals(abondance7_zi)
plot(sim_res_zi)

AIC(abondance7, abondance7_zi)
#pas mieux

# Tester si le problème vient des valeurs extrêmes
testOutliers(sim_res)

# Tester l'autocorrélation spatiale

# Agréger les résidus par lac
sim_res_agg <- recalculateResiduals(sim_res, group = data_model$Lake)

# Coordonnées moyennes par lac
coords <- data_model %>%
  group_by(Lake) %>%
  summarise(x = mean(Longitude), y = mean(Latitude))

# Test autocorrélation spatiale
testSpatialAutocorrelation(sim_res_agg, 
                           x = coords$x, 
                           y = coords$y)
#Le modèle est globalement valide
# La seule limitation est la légère déviation des résidus détectée par DHARMa 
# attribuable à la variabilité résiduelle entre lacs 
# ce que la forte variance de l'effet aléatoire Lake confirme

# Et le modèle de prévalence ? 

sim_res_prev <- simulateResiduals(prevalence5)
plot(sim_res_prev)
# meilleur que l'abondance au QQ plot

sim_res_prev_agg <- recalculateResiduals(sim_res_prev, group = data_model$Lake)
testSpatialAutocorrelation(sim_res_prev_agg, x = coords$x, y = coords$y)
# pas d'autocorrélation spatiale non plus

# Effet des grands cyclopoïdes sur l'abondance

pred <- ggpredict(abondance7, terms = "Cyclo_grand1_logstd")

ggplot(pred, aes(x = x, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 200)) +  # limite l'axe Y aux valeurs réalistes
  labs(x = "Grands cyclopoïdes (> 1mm) [standardisé]",
       y = "Abondance prédite d'Apophallus") +
  theme_classic()

data_lac <- data_model %>%
  group_by(Lake) %>%
  summarise(mean_para = mean(Parasites),
            mean_cyclo_grand1 = mean(Cyclo_grand1, na.rm = TRUE),
            .groups = "drop")

ggplot(data_lac, aes(x = mean_cyclo_grand1, y = mean_para)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  geom_label_repel(aes(label = Lake), size = 2.5) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           size = 4) +
  coord_cartesian(ylim = c(0, 140)) +
  labs(x = "Abondance moyenne des grands cyclopoïdes (> 1mm)",
       y = "Abondance moyenne d'Apophallus par lac") +
  theme_classic()

# Effet du Type de communauté sur l'abondance

data_model <- data_model %>%
  dplyr::mutate(Type_label = recode(Type, 
                             "Allo" = "Allopatrique",
                             "Symp" = "Sympatrique"))

ggplot(data_model, aes(x = Type_label, y = Parasites, fill = Type_label)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2, width = 0.2, size = 0.8) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_fill_manual(values = c("Allopatrique" = "#7BAFD4", 
                               "Sympatrique" = "#E8A87C")) +
  labs(x = "Type de communauté",
       y = "Abondance d'Apophallus") +
  theme_classic() +
  theme(legend.position = "none")

# Carte
data_model <- data_model %>%
  mutate(Longitude = -Longitude)

# Recréer data_carte
data_carte <- data_model %>%
  group_by(Lake, Latitude, Longitude) %>%
  summarise(mean_para = mean(Parasites),
            prevalence = mean(Presence),
            .groups = "drop")

# Télécharger le fond de carte
zone <- terra::ext(-74, -72.5, 46.4, 47.3)
fond <- get_tiles(zone, provider = "OpenStreetMap", zoom = 10)

# Graphique
carte_principale <- ggplot() +
  geom_spatraster_rgb(data = fond) +
  geom_point(data = data_carte, 
             aes(x = Longitude, y = Latitude, 
                 size = mean_para, color = mean_para),
             alpha = 0.8) +
  geom_label_repel(data = data_carte,
                   aes(x = Longitude, y = Latitude, label = Lake),
                   size = 2.5, max.overlaps = 20,
                   box.padding = 0.5) +
  scale_color_gradient(low = "#3182bd", high = "#de2d26",
                       name = "Abondance\nmoyenne") +
  scale_size_continuous(range = c(3, 10),
                        name = "Abondance\nmoyenne") +
  guides(color = guide_legend(), size = guide_legend()) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()

# Encart Québec
encart <- ggplot() +
  annotation_borders("world", regions = "Canada", 
                     fill = "lightgrey", colour = "white") +
  annotation_borders("state", fill = NA, colour = "white") +
  geom_rect(aes(xmin = -74, xmax = -72.5, ymin = 46.4, ymax = 47.3),
            fill = NA, color = "red", linewidth = 1) +
  coord_cartesian(xlim = c(-80, -60), ylim = c(44, 52)) +
  theme_void()

# Assembler
ggdraw(carte_principale) +
  draw_plot(encart, x = 0.25, y = 0.68, width = 0.22, height = 0.22)

# Taille cyclopoides par lacs
zoo %>%
  filter(object_annotation_category == "Cyclopoida") %>%
  ggplot(aes(x = reorder(Lake, Major_mm, median), y = Major_mm)) +
  geom_boxplot(fill = "#7BAFD4", outlier.shape = NA) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0.8, color = "orange", linetype = "dashed", linewidth = 0.8) +
  coord_cartesian(ylim = c(0, 1.6)) +
  labs(x = "Lac",
       y = "Taille des cyclopoïdes (mm)",
       caption = "Rouge = seuil 1mm | Orange = seuil 0.8mm") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculer la taille moyenne des cyclopoïdes par lac
taille_cyclo_lac <- zoo %>%
  filter(object_annotation_category == "Cyclopoida") %>%
  group_by(Lake) %>%
  summarise(mean_taille = mean(Major_mm, na.rm = TRUE),
            .groups = "drop")

# Joindre avec les parasites
data_taille <- left_join(data_carte, taille_cyclo_lac, by = "Lake")

# Graphique
ggplot(data_taille, aes(x = mean_taille, y = mean_para)) +
  geom_point(size = 3) +
  coord_cartesian(ylim = c(0, 140)) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  geom_label_repel(aes(label = Lake), size = 2.5) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  labs(x = "Taille moyenne des cyclopoïdes (mm)",
       y = "Abondance moyenne d'Apophallus") +
  theme_classic()

# Effet biomasse poissons sur cyclopoides ? 

data_model <- data_model %>%
  mutate(Biomasse_moy = as.numeric(gsub(",", ".", Biomasse_moy)))

data_lac3 <- data_model %>%
  group_by(Lake) %>%
  summarise(Cyclo_total = mean(Cyclo_total, na.rm = TRUE),
            Biomasse = mean(Biomasse_moy, na.rm = TRUE),
            .groups = "drop")

ggplot(data_lac3, aes(x = Biomasse, y = Cyclo_total)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "black") +
  geom_label_repel(aes(label = Lake), size = 2.5) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 4000, label.y = 450, size = 4) +
  labs(x = "Biomasse moyenne des poissons (g/trap)",
       y = "Abondance totale des cyclopoïdes") +
  theme_classic()

# Et avec les calanoides ? 

data_lac4 <- data_model %>%
  group_by(Lake) %>%
  summarise(Cala_total = mean(Cala_total, na.rm = TRUE),
            Biomasse = mean(Biomasse_moy, na.rm = TRUE),
            .groups = "drop")

ggplot(data_lac4, aes(x = Biomasse, y = Cala_total)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "black") +
  geom_label_repel(aes(label = Lake), size = 2.5) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 4000, label.y = max(data_lac4$Cala_total)*0.9, size = 4) +
  labs(x = "Biomasse moyenne des poissons (g/trap)",
       y = "Abondance totale des calanoides") +
  theme_classic()

# Effet de la superficie sur la prévalence

mean_sup <- mean(log(data_model$Superficie), na.rm = TRUE)
sd_sup <- sd(log(data_model$Superficie), na.rm = TRUE)

pred_sup$x_original <- exp(pred_sup$x * sd_sup + mean_sup)

ggplot(pred_sup, aes(x = x_original, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Superficie (ha)",
       y = "Probabilité d'infection par Apophallus") +
  theme_classic()

# Effet du Type
ggplot(pred_type, aes(x = x, y = predicted)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  scale_x_discrete(labels = c("Allo" = "Allopatrique", "Symp" = "Sympatrique")) +
  labs(x = "Type de communauté",
       y = "Probabilité d'infection par Apophallus") +
  theme_classic()

## Et après tout, si on regardait si les autres zooplanctons potentiellement
# carnivores peuvent être utiles à notre modèle ? 

zoo %>%
  filter(object_annotation_category %in% c("Chaoboridae", "Leptodora kindtii", "Acari", "Asplanchna")) %>%
  group_by(object_annotation_category) %>%
  summarise(n_total = n(),
            n_lacs = n_distinct(Lake),
            mean_par_lac = n_total / n_lacs)

# Créer les variables présence/absence par lac
# Pris un à un, cela va être compliqué d'avoir de la puissance statistique, 
# on va donc créer une autre variable "présence d'au moins un autre prédateur"

zoo_pa <- zoo %>%
  group_by(Lake) %>%
  summarise(
    Chaoborus_pa = as.integer(any(object_annotation_category == "Chaoboridae")),
    Asplanchna_pa = as.integer(any(object_annotation_category == "Asplanchna")),
    Acari_pa = as.integer(any(object_annotation_category == "Acari")),
    Leptodora_pa = as.integer(any(object_annotation_category == "Leptodora kindtii")),
    .groups = "drop"
  ) %>%
  mutate(Autres_predateurs = as.integer(Chaoborus_pa == 1 | 
                                          Asplanchna_pa == 1 | 
                                          Acari_pa == 1 | 
                                          Leptodora_pa == 1))

zoo_pa <- zoo_pa %>%
  filter(Lake %in% unique(data_model$Lake))

nrow(zoo_pa)
sum(zoo_pa$Autres_predateurs)

sum(zoo_pa$Autres_predateurs)
table(zoo_pa$Autres_predateurs)

data_model <- left_join(data_model, 
                        zoo_pa %>% select(Lake, Chaoborus_pa), 
                        by = "Lake")

# On prend finalement que Chaoborus parce que sinon on a une présence d'au moins
# un des 4 dans tous les lacs donc ça n'a pas de sens

# Modèle présence Chaoborus
abondance_chaoborus <- glmer.nb(Parasites ~ Cyclo_grand1_logstd + Chaoborus_pa + LF_std + (1|Lake),
                                data = data_model)

summary(abondance_chaoborus)
AIC(abondance7, abondance_chaoborus)

# Modèle abondance Asplanchna
asplanchna_ab <- zoo %>%
  group_by(Lake, ID) %>%
  summarise(
    Asplanchna_ab = sum(object_annotation_category == "Asplanchna"),
    .groups = "drop"
  ) %>%
  group_by(Lake) %>%
  summarise(
    Asplanchna_ab = mean(Asplanchna_ab),
    .groups = "drop"
  ) %>%
  filter(Lake %in% unique(data_model$Lake))

# Joindre
data_model <- left_join(data_model, asplanchna_ab, by = "Lake")

# Vérifier les NA
sum(is.na(data_model$Asplanchna_ab))

# Standardiser
data_model <- data_model %>%
  mutate(Asplanchna_logstd = as.numeric(scale(log(Asplanchna_ab + 1))))

# Modèle
abondance_asplanchna <- glmer.nb(Parasites ~ Asplanchna_logstd + Cyclo_grand1_logstd + LF_std + (1|Lake),
                                 data = data_model)

summary(abondance_asplanchna)
AIC(abondance7, abondance_asplanchna)
# Effet négatif donc dans la bonne direction biologiquement mais non significatif
# AIC plus élevée

# Pour abondance absolue Chao
chaoborus_ab <- zoo %>%
  group_by(Lake, ID) %>%
  summarise(
    Chaoborus_ab = sum(object_annotation_category == "Chaoboridae"),
    .groups = "drop"
  ) %>%
  group_by(Lake) %>%
  summarise(
    Chaoborus_ab = mean(Chaoborus_ab),
    .groups = "drop"
  ) %>%
  filter(Lake %in% unique(data_model$Lake))

# Joindre et standardiser
data_model <- left_join(data_model, chaoborus_ab, by = "Lake") %>%
  mutate(Chaoborus_logstd = as.numeric(scale(log(Chaoborus_ab + 1))))

# Modèle
abondance_chaoborus <- glmer.nb(Parasites ~ Chaoborus_logstd + Cyclo_grand1_logstd + LF_std + (1|Lake),
                                data = data_model)

summary(abondance_chaoborus)
AIC(abondance7, abondance_chaoborus)
