
ordered <- function(unordered, ordered_by, logical){
      unordered[order(unordered[ordered_by], decreasing = logical),]
}

ordered.dem <- ordered(demographics, "age", FALSE)


#summarize all data and look for strange things/figure out what was recorded for different types of death
summary(demographics)
summary(demographics$cause_of_death)

#885 animals have NA listed for "was_control_mock_treated" and all are JM11
summary(demographics$was_control_mock_treated)
demographics[is.na(demographics$was_control_mock_treated),]

#42 animals were first irradiated after their age at death... makes no sense
summary(demographics$first_irrad > demographics$age)
demographics[demographics$first_irrad > demographics$age,]
summary(demographics[demographics$first_irrad > demographics$age,])
#output these animals to get list of animal IDs and look up more info on mice
list.weird.animals <- demographics[demographics$first_irrad > demographics$age,]
write.csv(list.weird.animals, file = "death_before_irradiation.csv",
          row.names=TRUE)

#age is listed for all animals
summary(is.na(demographics$age))

# Accidental death info / mix of conditions
summary(subset(demographics, cause_of_death == "Accidental death"))

# Escaped during irradiation info / mix of conditions
summary(subset(demographics, cause_of_death == "Escaped during irradiation"))

# Discard info / mix of conditions / 150 NA for mock treated... seems high
summary(subset(demographics, cause_of_death == "Discard"))

# Grahn mice, breeder / all expt 13, all 60 fractions, all died around 550 days / mix of dose and type of radn
summary(subset(demographics, cause_of_death == "Grahn mice, breeder"))

# Improper irradiation info / mix of conditions
summary(subset(demographics, cause_of_death == "Improper irradiation"))

# Missing info / mix of conditions
summary(subset(demographics, cause_of_death == "Missing"))

# Removal to another experiment info / mix of conditions
summary(subset(demographics, cause_of_death == "Removal to another experiment"))

# Sacrifice, programmed info / only experiment 2 and 4, mix of conditions otherwise
summary(subset(demographics, cause_of_death == "Sacrifice, programmed"))

###################################### data I want to use ########################################
# Sacrifice, moribund info / mix of conditions 
summary(subset(demographics, cause_of_death == "Sacrifice, moribund"))
head(subset(demographics, cause_of_death == "Sacrifice, moribund"))

# Died info / mix of conditions 
summary(subset(demographics, cause_of_death == "Died"))


######## Check for mice that don't have lethal cause of death listed #######
lethal_columns <- grep("_L", names(all_info))
all_info$lethal_sum <- rowSums(all_info[,lethal_columns])
summary(all_info$lethal_sum == 0)

mice_to_use <- mice_to_use[mice_to_use$lethal_sum == 1,] 

#psuedo code
#combine demographics and macros for all
#find lethal columbs and sum it
#find out if any are missing one