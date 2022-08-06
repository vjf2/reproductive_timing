get_age <- function(data,
                    column = "Dolphin.ID",
                    atDate = "Conception.Date",
                    lhdata) {

  data[, "birth_date"] <- lhdata$Birth.Date[match(data[, column], lhdata[, "Dolphin.ID"])]
  data[, "birth_date"] <- as.Date(data[, "birth_date"])

  age <- data[, atDate] - data[, "birth_date"]
  age <- as.numeric(age) / 365.25
  return(age)

}
