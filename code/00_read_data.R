rm(list=ls())

library(tidyverse)

person <- read.csv("data/raw/FARS2024NationalCSV/person.csv")
str(person)

vehicle <- read.csv("data/raw/FARS2024NationalCSV/vehicle.csv")
str(vehicle)

accident <- read.csv("data/raw/FARS2024NationalCSV/accident.csv")
str(accident)

unique(accident$YEAR)
