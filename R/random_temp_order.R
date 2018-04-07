setwd("~/Dropbox/Ldeli_inc_metab/")

data <- data.frame(batch = rep(c(1,2), each = 4),
                   day = rep(c(1:4), 2))
                   
data$temp <- unlist(list((replicate(2,sample(c(20,24,28,32), replace = F)))))

write.csv(data, "data/ld_temporder.csv", row.names = F)
