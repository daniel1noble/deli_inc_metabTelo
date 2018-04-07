#Script for randomly selection 25 hot and 25 cold babies for the metabolic rate and telomere study from the 2017 -2018 F1 cohort

setwd("~/Dropbox/Ldeli_inc_metab")

library(dplyr)
library(ggplot2)

dat <- read.csv("data/raw/season2_hatchlings.csv")
dat$egg_incub_temp <- as.factor(dat$egg_incub_temp)
dat$chrono_age <- as.Date(as.character(Sys.Date())) - as.Date(as.character(dat$liz_hatch_date), format = "%d/%m/%Y")
dat$chrono_age <- as.numeric(dat$chrono_age)

str(dat)
names(dat)[1] <- "liz_id"
summary(dat)

ggplot(dat, aes(chrono_age)) + 
  geom_histogram(alpha = 0.2, binwidth = 8) +
  theme_bw()

#Sampling from our two treatments
set.seed(106)

cold_samp <- sample(dat[dat$egg_incub_temp == "23", "liz_id"], 25)
hot_samp <- sample(dat[dat$egg_incub_temp == "29", "liz_id"], 25)

my_treat_ids <- c(as.vector(cold_samp) , as.vector(hot_samp)) 

final_dat <- rbind(dat[dat$liz_id %in% my_treat_ids,])

final_dat %>% filter(egg_incub_temp == "23") %>% select(liz_id) %>% unique
final_dat %>% filter(egg_incub_temp == "29") %>% select(liz_id) %>% unique

#age distribution

ggplot(final_dat, aes(y = chrono_age, x = factor(egg_incub_temp))) + 
  geom_boxplot(alpha = 0.2)

ggplot(final_dat, aes(chrono_age, fill = factor(egg_incub_temp))) + 
  geom_histogram(alpha = 0.2, binwidth = 8) + 
  scale_fill_manual(values=c("#4DAF4A", "#984EA3")) + 
  theme_bw()

#Creating randomisation dataframe for 10 sampling periods, two batches, 2 incubators and 25 chambers spread across both incbuators

randat <- data.frame(period = rep(c(1:10), each = 50),
                     batch = rep(c(1:2), times = c(25, 25)),
                     incubator = rep(c(1,2,1,2), times = c(12,13,12,13)),
                     chamber = rep(c(1:25, 1:25), times = c(10)))

randat <- randat[order(randat$batch),]    

#Assigning babies to batches
#Assign to two difference batches
#This confirms that batch 1 and batch 2 ids are not the same. Good
batch1_ids <- sample(my_treat_ids, size = 25, replace = F) %>% sort
table(final_dat[final_dat$liz_id %in% batch1_ids, "egg_incub_temp"])

[1] "ld0449" "ld0451" "ld0455" "ld0457" "ld0460" "ld0461" "ld0471" "ld0472"
[9] "ld0473" "ld0479" "ld0481" "ld0484" #= Dead "ld0494" "ld0497" "ld0503" "ld0507"
[17] "ld0512" "ld0513" "ld0519" "ld0539" "ld0543" "ld0544" "ld0546" "ld0550"
[25] "ld0553"

batch2_ids <- my_treat_ids[! my_treat_ids %in% batch1_ids] %>% sort 
table(final_dat[final_dat$liz_id %in% batch2_ids, "egg_incub_temp"])

[1] "ld0458" "ld0463" "ld0466" "ld0469" "ld0477" "ld0487" "ld0488" "ld0489"
[9] "ld0490" "ld0492" "ld0493" "ld0502" "ld0504" "ld0505" "ld0510" "ld0511"
[17] "ld0515" "ld0521" "ld0527" "ld0528" "ld0529" "ld0530" "ld0540" "ld0541"
[25] "ld0545"

#Now I want my lizards to remain in the same batch and the order can be mixed so each baby can be in a diff chamber every session
#Sampling B1 ids for 10 periods and randomise them into different order 10 times so they are spread across chambers
allb1 <- replicate(10,sample(batch1_ids, replace = F))
allb1 <- unlist(allb1)
length(allb1)

table(final_dat[final_dat$liz_id %in% unique(unlist(allb1)), "egg_incub_temp"])
table(final_dat[final_dat$liz_id %in% unique(batch1_ids), "egg_incub_temp"])

allb2 <- replicate(10,sample(batch2_ids, replace = F))
allb2 <- unlist(allb2)
length(allb2)

summary(left_join(data.frame(liz_id = unique(unlist(batch2_ids))), final_dat[,1:2]))
summary(left_join(data.frame(liz_id = batch2_ids), final_dat[,1:2]))

table(final_dat[final_dat$liz_id %in% unique(unlist(allb2)), "egg_incub_temp"])
table(final_dat[final_dat$liz_id %in% unique(batch2_ids), "egg_incub_temp"])

randat$liz_id <- c(allb1, allb2)
randat <- randat[order(randat$batch),]    

#Merge with lizard details

final_dat2 <-left_join(randat, final_dat[,1:2])
write.csv(final_dat2, row.names = F, "data/final/Tinc_hatchlings.csv")

#I want unique ids with egg treatments for both batches
final_dat2 %>% filter(batch == "1") %>% select(liz_id, egg_incub_temp) %>% unique %>%summary
final_dat2 %>% filter(batch == "2") %>% select(liz_id, egg_incub_temp) %>% unique %>%summary

#This confirms that batch 1 lizards are always batch 1 lizards and same for batch 2
final_dat2 %>% filter(batch == "1" & period == "1") %>% pull(liz_id) %>% sort == final_dat2 %>% filter(batch == "1" & period == "2") %>% pull(liz_id) %>% sort
final_dat2 %>% filter(batch == "2" & period == "1") %>% pull(liz_id) %>% sort == final_dat2 %>% filter(batch == "2" & period == "2") %>% pull(liz_id) %>% sort

#This confirms that the order of lizards across chambers are different which is what I want
final_dat2 %>% filter(batch == "1" & period == "5") %>% select(liz_id, chamber) 
final_dat2 %>% filter(batch == "1" & period == "10") %>% select(liz_id, chamber) 

final_dat2 %>% filter(batch == "2" & period == "5") %>% select(liz_id, chamber) 
final_dat2 %>% filter(batch == "2" & period == "10") %>% select(liz_id, chamber) 

#Descriptive stats of lizards in batch and incubator is it balanced? 
final_dat2 %>% filter(batch == "1") %>% select(liz_id, egg_incub_temp) %>% unique %>% summary
final_dat2 %>% filter(batch == "2") %>% select(liz_id, egg_incub_temp) %>% unique %>% summary   

final_dat2 %>% filter(batch == "1" & incubator == "1" & period == "1") %>% select(liz_id, egg_incub_temp) %>% summary()
final_dat2 %>% filter(batch == "1" & incubator == "2" & period == "1") %>% select(liz_id, egg_incub_temp) %>% summary()
final_dat2 %>% filter(batch == "1" & period == "1") %>% select(liz_id, egg_incub_temp) %>% summary()


final_dat2 %>% filter(batch == "2" & incubator == "1" & period == "1") %>% select(liz_id, egg_incub_temp) %>% summary()
final_dat2 %>% filter(batch == "2" & incubator == "2" & period == "1") %>% select(liz_id, egg_incub_temp) %>% summary()
final_dat2 %>% filter(batch == "2" & period == "1") %>% select(liz_id, egg_incub_temp) %>% summary()

##Constructing randomisation dataframe for temperatures to incbuators
#new way after chatting to Dan

my.temp.dat <- data.frame(batch = rep(c(1,2), each = 40),
                          period = rep(c(1:10), each= 2),
                          day_id = rep(c(1,2), 40, each = 1),
                          incb_id = rep(c(1,2), each = 20),
                          temp_order = unlist(list((replicate(40, sample(c(0,1)))))),
                          t1 = rep("NA", 80),
                          t2 = rep("NA", 80))

low_temps <- c(24, 26, 28)
high_temps <- c(30, 32, 34)

my.temp.dat$t1 <- as.numeric(my.temp.dat$t1)
my.temp.dat$t2 <- as.numeric(my.temp.dat$t2)

#Have to do this by hand I think
head(my.temp.dat)
View(my.temp.dat)

#batch 1, incubator 1

#p1 d1
my.temp.dat[1,"t1"] <- sample(high_temps, 1)
my.temp.dat[1,"t2"] <- sample(low_temps, 1)
#p1 d2
my.temp.dat[2,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[1,"t2"]], 1)
my.temp.dat[2,"t2"] <-  sample(high_temps[!high_temps %in% my.temp.dat[1,"t1"]], 1)

#b1, p2, d1, i1 
my.temp.dat[3,"t1"] <- sample(high_temps, 1)
my.temp.dat[3,"t2"] <- sample(low_temps, 1)
#b1, p2, d2, i1
my.temp.dat[4,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[3,"t2"]], 1)
my.temp.dat[4,"t2"] <-  sample(high_temps[!high_temps %in% my.temp.dat[3,"t1"]], 1)

#b1, p3, d1, i1 
my.temp.dat[5,"t1"] <- sample(high_temps, 1)
my.temp.dat[5,"t2"] <- sample(low_temps, 1)
#b1, p3, d2, i1
my.temp.dat[6,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[5,"t2"]], 1)
my.temp.dat[6,"t2"] <-  sample(high_temps[!high_temps %in% my.temp.dat[5,"t1"]], 1)

#b1, p4, d1, i1 
my.temp.dat[7,"t1"] <- sample(low_temps, 1)
my.temp.dat[7,"t2"] <- sample(high_temps, 1)
#b1, p4, d2, i1
my.temp.dat[8,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[7,"t2"]], 1)
my.temp.dat[8,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[7,"t1"]], 1)

#b1, p5, d1, i1 
my.temp.dat[9,"t1"] <- sample(low_temps, 1)
my.temp.dat[9,"t2"] <- sample(high_temps, 1)
#b1, p5, d2, i1
my.temp.dat[10,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[9,"t2"]], 1)
my.temp.dat[10,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[9,"t1"]], 1)

#b1, p6, d1, i1 
my.temp.dat[11,"t1"] <- sample(high_temps, 1)
my.temp.dat[11,"t2"] <- sample(low_temps, 1)
#b1, p6, d2, i1
my.temp.dat[12,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[11,"t2"]], 1)
my.temp.dat[12,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[11,"t1"]], 1)

#b1, p7, d1, i1 
my.temp.dat[13,"t1"] <- sample(high_temps, 1)
my.temp.dat[13,"t2"] <- sample(low_temps, 1)
#b1, p7, d2, i1
my.temp.dat[14,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[13,"t2"]], 1)
my.temp.dat[14,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[13,"t1"]], 1)

#b1, p8, d1, i1 
my.temp.dat[15,"t1"] <- sample(low_temps, 1)
my.temp.dat[15,"t2"] <- sample(high_temps, 1)
#b1, p8, d2, i1
my.temp.dat[16,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[15,"t2"]], 1)
my.temp.dat[16,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[15,"t1"]], 1)

#b1, p9, d1, i1 
my.temp.dat[17,"t1"] <- sample(low_temps, 1)
my.temp.dat[17,"t2"] <- sample(high_temps, 1)
#b1, p9, d2, i1
my.temp.dat[18,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[17,"t2"]], 1)
my.temp.dat[18,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[17,"t1"]], 1)

#b1, p10, d1, i1 
my.temp.dat[19,"t1"] <- sample(low_temps, 1)
my.temp.dat[19,"t2"] <- sample(high_temps, 1)
#b1, p10, d2, i1
my.temp.dat[20,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[19,"t2"]], 1)
my.temp.dat[20,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[19,"t1"]], 1)

#batch 1, incubator 2
#p1, d1
my.temp.dat[21,"t1"] <- sample(low_temps, 1)
my.temp.dat[21,"t2"] <- sample(high_temps, 1)
#p1, d2
my.temp.dat[22,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[21,"t2"]], 1)
my.temp.dat[22,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[21,"t1"]], 1)
#p2, d1
my.temp.dat[23,"t1"] <- sample(low_temps, 1)
my.temp.dat[23,"t2"] <- sample(high_temps, 1)
#p2, d2
my.temp.dat[24,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[23,"t2"]], 1)
my.temp.dat[24,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[23,"t1"]], 1)
#p3, d1
my.temp.dat[25,"t1"] <- sample(low_temps, 1)
my.temp.dat[25,"t2"] <- sample(high_temps, 1)
#p3, d2
my.temp.dat[26,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[25,"t2"]], 1)
my.temp.dat[26,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[25,"t1"]], 1)
#p4, d1
my.temp.dat[27,"t1"] <- sample(high_temps, 1)
my.temp.dat[27,"t2"] <- sample(low_temps, 1)
#p4, d2
my.temp.dat[28,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[27,"t2"]], 1)
my.temp.dat[28,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[27,"t1"]], 1)
#p5, d1
my.temp.dat[29,"t1"] <- sample(low_temps, 1)
my.temp.dat[29,"t2"] <- sample(high_temps, 1)
#p5, d2
my.temp.dat[30,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[29,"t2"]], 1)
my.temp.dat[30,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[29,"t1"]], 1)
#p6, d1
my.temp.dat[31,"t1"] <- sample(low_temps, 1)
my.temp.dat[31,"t2"] <- sample(high_temps, 1)
#p6, d2
my.temp.dat[32,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[31,"t2"]], 1)
my.temp.dat[32,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[31,"t1"]], 1)
#p7, d1
my.temp.dat[33,"t1"] <- sample(low_temps, 1)
my.temp.dat[33,"t2"] <- sample(high_temps, 1)
#p7, d2
my.temp.dat[34,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[33,"t2"]], 1)
my.temp.dat[34,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[33,"t1"]], 1)
#p8, d1
my.temp.dat[35,"t1"] <- sample(low_temps, 1)
my.temp.dat[35,"t2"] <- sample(high_temps, 1)
#p8, d2
my.temp.dat[36,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[35,"t2"]], 1)
my.temp.dat[36,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[35,"t1"]], 1)
#p9, d1
my.temp.dat[37,"t1"] <- sample(high_temps, 1)
my.temp.dat[37,"t2"] <- sample(low_temps, 1)
#p9, d2
my.temp.dat[38,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[37,"t2"]], 1)
my.temp.dat[38,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[37,"t1"]], 1)
#p10, d1
my.temp.dat[39,"t1"] <- sample(high_temps, 1)
my.temp.dat[39,"t2"] <- sample(low_temps, 1)
#p10, d2
my.temp.dat[40,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[39,"t2"]], 1)
my.temp.dat[40,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[39,"t1"]], 1)

#batch 2 incubator 1
#p1 day 1
my.temp.dat[41,"t1"] <- sample(high_temps, 1)
my.temp.dat[41,"t2"] <- sample(low_temps, 1)
#p1, d2
my.temp.dat[42,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[41,"t2"]], 1)
my.temp.dat[42,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[41,"t1"]], 1)
#p2 day 1
my.temp.dat[43,"t1"] <- sample(high_temps, 1)
my.temp.dat[43,"t2"] <- sample(low_temps, 1)
#p2, d2
my.temp.dat[44,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[43,"t2"]], 1)
my.temp.dat[44,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[43,"t1"]], 1)
#p3, d1
my.temp.dat[45,"t1"] <- sample(low_temps, 1)
my.temp.dat[45,"t2"] <- sample(high_temps, 1)
#p3, d2
my.temp.dat[46,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[45,"t2"]], 1)
my.temp.dat[46,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[45,"t1"]], 1)
#p4, d1
my.temp.dat[47,"t1"] <- sample(low_temps, 1)
my.temp.dat[47,"t2"] <- sample(high_temps, 1)
#p4, d2
my.temp.dat[48,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[47,"t2"]], 1)
my.temp.dat[48,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[47,"t1"]], 1)
#p5, d1
my.temp.dat[49,"t1"] <- sample(low_temps, 1)
my.temp.dat[49,"t2"] <- sample(high_temps, 1)
#p5, d2
my.temp.dat[50,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[49,"t2"]], 1)
my.temp.dat[50,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[49,"t1"]], 1)
#p6 day 1
my.temp.dat[51,"t1"] <- sample(high_temps, 1)
my.temp.dat[51,"t2"] <- sample(low_temps, 1)
#p6, d2
my.temp.dat[52,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[51,"t2"]], 1)
my.temp.dat[52,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[51,"t1"]], 1)
#p7, d1
my.temp.dat[53,"t1"] <- sample(low_temps, 1)
my.temp.dat[53,"t2"] <- sample(high_temps, 1)
#p7, d2
my.temp.dat[54,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[53,"t2"]], 1)
my.temp.dat[54,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[53,"t1"]], 1)
#p8 day 1
my.temp.dat[55,"t1"] <- sample(high_temps, 1)
my.temp.dat[55,"t2"] <- sample(low_temps, 1)
#p8, d2
my.temp.dat[56,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[55,"t2"]], 1)
my.temp.dat[56,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[55,"t1"]], 1)
#p9 day 1
my.temp.dat[57,"t1"] <- sample(high_temps, 1)
my.temp.dat[57,"t2"] <- sample(low_temps, 1)
#p9, d2
my.temp.dat[58,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[57,"t2"]], 1)
my.temp.dat[58,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[57,"t1"]], 1)
#p10, d1
my.temp.dat[59,"t1"] <- sample(low_temps, 1)
my.temp.dat[59,"t2"] <- sample(high_temps, 1)
#p10, d2
my.temp.dat[60,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[59,"t2"]], 1)
my.temp.dat[60,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[59,"t1"]], 1)

#batch 2 incubator 2
#p1, d1
my.temp.dat[61,"t1"] <- sample(low_temps, 1)
my.temp.dat[61,"t2"] <- sample(high_temps, 1)
#p1, d2
my.temp.dat[62,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[61,"t2"]], 1)
my.temp.dat[62,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[61,"t1"]], 1)
#p2, d1
my.temp.dat[63,"t1"] <- sample(low_temps, 1)
my.temp.dat[63,"t2"] <- sample(high_temps, 1)
#p2, d2
my.temp.dat[64,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[63,"t2"]], 1)
my.temp.dat[64,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[63,"t1"]], 1)
#p3, d1
my.temp.dat[65,"t1"] <- sample(low_temps, 1)
my.temp.dat[65,"t2"] <- sample(high_temps, 1)
#p3, d2
my.temp.dat[66,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[65, "t2"]], 1)
my.temp.dat[66,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[65,"t1"]], 1)
#p4, d1
my.temp.dat[67,"t1"] <- sample(low_temps, 1)
my.temp.dat[67,"t2"] <- sample(high_temps, 1)
#p4, d2
my.temp.dat[68,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[67,"t2"]], 1)
my.temp.dat[68,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[67,"t1"]], 1)
#p5, d1
my.temp.dat[69,"t1"] <- sample(low_temps, 1)
my.temp.dat[69,"t2"] <- sample(high_temps, 1)
#p5, d2
my.temp.dat[70,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[69,"t2"]], 1)
my.temp.dat[70,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[69,"t1"]], 1)
#p6, d1
my.temp.dat[71,"t1"] <- sample(low_temps, 1)
my.temp.dat[71,"t2"] <- sample(high_temps, 1)
#p6, d2
my.temp.dat[72,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[71,"t2"]], 1)
my.temp.dat[72,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[71,"t1"]], 1)
#p7, d1
my.temp.dat[73,"t1"] <- sample(low_temps, 1)
my.temp.dat[73,"t2"] <- sample(high_temps, 1)
#p7, d2
my.temp.dat[74,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[73,"t2"]], 1)
my.temp.dat[74,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[73,"t1"]], 1)
#p8, d1
my.temp.dat[75,"t1"] <- sample(low_temps, 1)
my.temp.dat[75, "t2"] <- sample(high_temps, 1)
#p8, d2
my.temp.dat[76,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[75,"t2"]], 1)
my.temp.dat[76,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[75,"t1"]], 1)
#p9, d1
my.temp.dat[77,"t1"] <- sample(low_temps, 1)
my.temp.dat[77,"t2"] <- sample(high_temps, 1)
#p9, d2
my.temp.dat[78,"t1"] <- sample(high_temps[!high_temps %in% my.temp.dat[77,"t2"]], 1)
my.temp.dat[78,"t2"] <- sample(low_temps[!low_temps %in% my.temp.dat[77,"t1"]], 1)
#p10, day 1
my.temp.dat[79,"t1"] <- sample(high_temps, 1)
my.temp.dat[79,"t2"] <- sample(low_temps, 1)
#p10, d2
my.temp.dat[80,"t1"] <- sample(low_temps[!low_temps %in% my.temp.dat[79,"t2"]], 1)
my.temp.dat[80,"t2"] <- sample(high_temps[!high_temps %in% my.temp.dat[79,"t1"]], 1)

write.csv(my.temp.dat, row.names = F, "data/final/Tinc_randtemp_V2.csv")

#### Archived code below
my.temp.dat[my.temp.dat$temp_order == "0" & my.temp.dat$day_id == "1", "t1"] <- as.vector(unlist(list(replicate(40,sample(low_temps, 1, replace  = F)))))

my.temp.dat[my.temp.dat$temp_order == "0" & my.temp.dat$day_id == "1","t2"] <- as.vector(unlist(list(replicate(40,sample(high_temps, 1, replace  = F)))))

View(my.temp.dat)

#old way
my.temps <- c(seq(24,34, 2))
my.temp.dat <- data.frame(batch = rep(c(1,2), each = 80),
                          period = rep(c(1:10), each = 8),
                          day_id = rep(c(1,2), 20, each = 2),
                          temp_id = rep(c(1:4), 20),
                          incb_id = rep(c(1,2), each = 4),
                          temp = unlist(list(replicate(40,sample(my.temps, 4, replace = F)))))
table(my.temp.dat$temp) #8 measures at a temp in total should be about ~13 or 14 measures per temp

write.csv(my.temp.dat, row.names = F, "data/randtemp.csv")

my.temp.dat %>% filter(batch == "1" , incb_id == "1") %>% pull(temp) %>% table
#need to shift 
#1 x 28 to 32 #period 3
#1 x 30 to 24 # period 4 

my.temp.dat %>% filter(batch == "1" , incb_id == "2") %>% pull(temp) %>% table
#need to shift 
#3 x 28 to 30 #period 2, 4, 8 
#1 x 34 to 26 #period 2

my.temp.dat %>% filter(batch == "2", incb_id == "1") %>% pull(temp) %>% table #all good 
#1 x 24 to 34 #period 4

my.temp.dat %>% filter(batch == "2", incb_id == "2") %>% pull(temp) %>% table #all good 
#1 x 30 to 24 #period 1
#1 x 26 to 30 #period 5

my.temp.dat <- read.csv("data/randtemp.csv") #read it back in to check, looks good


#test lizards I can use
#Will not contain my hot or cold_samp ids

set.seed(567)

test_ids <- sample(dat[!dat$liz_id %in% my_treat_ids,"liz_id"], 25)
sort(test_ids)

test_dat <- dat[dat$liz_id %in% test_ids,]

ggplot(test_dat, aes(y = chrono_age, x = factor(egg_incub_temp))) + 
  geom_boxplot(alpha = 0.2)

ggplot(test_dat, aes(chrono_age, fill = factor(egg_incub_temp))) + 
  geom_histogram(alpha = 0.2, binwidth = 8) + 
  scale_fill_manual(values=c("#4DAF4A", "#984EA3")) + 
  theme_bw()

summary(test_dat)

write.csv(test_dat, row.names = F, "data/Tinc_test_hatchlings.csv")
test_babies <- read.csv("data/Tinc_test_hatchlings.csv")
set.seed(757)

test_babies <- test_babies%>% filter(! liz_id == "ld0483")

brains <- rbind(test_babies %>% filter(egg_incub_temp == "23") %>% sample_n(5),
      test_babies %>% filter(egg_incub_temp == "29") %>% sample_n(5))

write.csv(brains, "~/Google Drive/2 - Projects/Ldeli_tinc_cogntion/brains_s2_f1s.csv")

brain_ids <- brains$liz_id
length(my_treat_ids)
length(brain_ids)

#Whose left for Tpref experiment? 
dat %>% filter(!liz_id %in% my_treat_ids) %>% filter(!liz_id %in% brain_ids)
dat %>% filter(!liz_id %in% my_treat_ids) %>% filter(!liz_id %in% brain_ids) %>% tally(egg_incub_temp == "23") #20 cold lizards
dat %>% filter(!liz_id %in% my_treat_ids) %>% filter(!liz_id %in% brain_ids) %>% tally(egg_incub_temp == "29") #14 cold lizards

write.csv(dat %>% filter(!liz_id %in% my_treat_ids) %>% filter(!liz_id %in% brain_ids) %>% as.data.frame(), row.names = F, "~/Google Drive/2 - Projects/Ldeli_tinc_cogntion/S2_F1_whoseleft_Tpref.csv")

dat %>% filter(!liz_id %in% my_treat_ids) %>% filter(!liz_id %in% brain_ids) %>% sample_n(4)
