### Code by: Bradley Rentz
### R Script for PaPE Conference poster: "And finally, no geminates!
### Pohnpeian consonantal length contrasts in initial, medial, and final
### position", Bradley Rentz and Victoria Anderson, 2017
### github.com/rentzb/pape
### author website: http://rentz.weebly.com
### Mac OS 10.12.4
### Base R version 3.3.3 "Another Canoe"

### dplyr version 0.5.0
### ggplot2 version 2.2.1
### rstanarm version 2.15.3
### yarrr version 0.1.5
### bayesplot version 1.2.0


library(ggplot2)
library(dplyr)
library(tidyr)
library(rstanarm)
library(plotly)
library(yarrr)
library(bayesplot)


options (mc.cores=parallel::detectCores ()) # Run on multiple cores

set.seed (3875)

setwd("~/Documents/UH/QP/qp1_data/results") # change this to your wd

### Pre-processing and cleaning
#load geminate data (303 obs)
geminates <- read.csv("geminates.csv")

# remove outliers
#remove 3 sd or greater for length grouped by speaker, consonant, and c_length (0 removed)
geminates.3sd <- geminates %>%
  group_by(speaker,consonant,word_place,c_length) %>%
  mutate(plusthreesdmean = mean(length)+3*sd(length),minusthreesd=mean(length)-3*sd(length))
geminates.clean <-geminates.3sd %>%
  filter(length < plusthreesdmean) %>%
  filter(length>minusthreesd)

# values over 0.35s are way too long so must be measurement error so remove
geminates.clean <- geminates.clean %>%
  filter(length< 0.35) # removes 3 out of 303 obs so 0.99%

geminates.clean$frame <- as.factor(geminates.clean$frame)

#data subset by place in word
geminates.clean.initial <- geminates.clean[geminates.clean$word_place=="initial",]
geminates.clean.medial <- geminates.clean[geminates.clean$word_place=="medial",]
geminates.clean.final <- geminates.clean[geminates.clean$word_place=="final",]

### Regression analyses by place in word
## note: depending on system, regression can be somewhat slow

# word initial

# note: multiply length*1000 to convert from seconds to ms
geminates.clean.initial.blmer = stan_lmer(length*1000 ~ consonant * c_length * frame + (1+consonant*c_length+frame|speaker) + (1+consonant*c_length+frame|word), data=geminates.clean.initial,
                                          prior_intercept = normal(0, 20),
                                          prior = normal(0, 40),
                                          prior_covariance = decov(regularization = 2),
                                          chains = 4,
                                          iter = 2000,adapt_delta=0.99999999999999999)

launch_shinystan(geminates.clean.initial.blmer) # use this to access shinystan diagnostic page

plot(geminates.clean.initial.blmer,pars=c("(Intercept)","c_lengthshort","consonantmw"))

sink("sink-summary-blmer-geminate-initial.txt")
print(summary(geminates.clean.initial.blmer),digits=2)
sink()

posterior_initial <- as.matrix(geminates.clean.initial.blmer)

plot_title <- ggtitle("BHLM for Word Initial Geminate Pairs","with 95% credible intervals")
mcmc_areas(posterior_initial,
           pars = c("(Intercept)","c_lengthshort","consonantmw","frame2","frame3"),
           prob=0.95) + plot_title + ggplot2::xlab("Duration (ms)")

# word medial

geminates.clean.medial.blmer = stan_lmer(length*1000 ~ consonant * c_length * frame + (1+consonant*c_length+frame|speaker) + (1+consonant*c_length+frame|word), data=geminates.clean.medial,
                                         prior_intercept = normal(0, 20),
                                         prior = normal(0, 40),
                                         prior_covariance = decov(regularization = 2),
                                         chains = 4,
                                         iter = 2000,adapt_delta=0.9999)

launch_shinystan(geminates.clean.medial.blmer) # use this to access shinystan diagnostic page


plot(geminates.clean.medial.blmer,pars=c("(Intercept)","c_lengthshort","consonantm","consonantmw","consonantn","consonantr"))

sink("sink-summary-blmer-geminate-medial.txt") # saves regression summary output as .txt file
print(summary(geminates.clean.medial.blmer),digits=2)
sink()

posterior_medial <- as.matrix(geminates.clean.medial.blmer)

plot_title <- ggtitle("BHLM for Word Medial Geminate Pairs","with 95% credible intervals")
mcmc_areas(posterior_medial,
           pars = c("(Intercept)","c_lengthshort","consonantm","consonantmw","consonantn","consonantr","frame2","frame3"),
           prob=0.95) + plot_title + ggplot2::xlab("Duration (ms)")


# word final

geminates.clean.final.blmer = stan_lmer(length*1000 ~ consonant * c_length * frame + (1+consonant*c_length+frame|speaker) + (1+consonant*c_length+frame|word), data=geminates.clean.final,
                                        prior_intercept = normal(0, 20),
                                        prior = normal(0, 40),
                                        prior_covariance = decov(regularization = 2),
                                        chains = 4,
                                        iter = 2000,adapt_delta=0.99999999999999)

launch_shinystan(geminates.clean.final.blmer) # use this to access shinystan diagnostic page


plot(geminates.clean.final.blmer,pars=c("(Intercept)","c_lengthshort","consonantmw"))

sink("sink-summary-blmer-geminate-final.txt")
print(summary(geminates.clean.final.blmer),digits=2)
sink()
posterior_final <- as.matrix(geminates.clean.final.blmer)

plot_title <- ggtitle("BHLM for Word Final Geminate Pairs","with 95% credible intervals")
mcmc_areas(posterior_final,
           pars = c("(Intercept)","c_lengthshort","consonantmw","frame2","frame3"),
           prob=0.95) + plot_title + ggplot2::xlab("Duration (ms)")



### plots

# change order of levels for word_place
geminates.clean$word_place=factor(geminates.clean$word_place,levels=c("initial","medial","final"))

# rdi plot (note: in some version of yarrr there is a bug causing some of the consonants to not be plotted. to fix, make separate version of plot for each place in word)
pirateplot(length*1000~c_length+consonant+word_place,data=geminates.clean,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="Duration (ms)",gl.col = "white",xlab="Consonant",cex.lab=0.4)


### Summary statistics
### means for geminate pairs

geminates.means <- geminates.clean %>%
  group_by(consonant, c_length,word_place) %>%
  summarise(cmean = mean(length), stderr = sd(length)/sqrt(length(length)))

geminates.meansoverall <- geminates.clean %>%
  group_by(consonant, c_length) %>%
  summarise(means = mean(length))

