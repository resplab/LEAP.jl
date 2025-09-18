# mortality rate

library(tidyverse)
library(here)
min_year <- 1996
max_year <- 2020
Canada_row <- c(22,111)

year_index <- function(year){
  year-1980+1+3
}

CA <- NULL
BC <- NULL

for(year in min_year:max_year){
  life_table_male <- openxlsx::read.xlsx(
      xlsxFile=here("data_generation/public_dataset/1980-2020_Tbl_1YR-eng.xlsx"),
      sheet=year_index(year),
      startRow=22) %>% 
    select(1,4,5)
  
  life_table_female <- openxlsx::read.xlsx(
      xlsxFile=here("data_generation/public_dataset/1980-2020_Tbl_1YR-eng.xlsx"),
      sheet=year_index(year),
      startRow=22) %>% 
    select(11,14,15)
  
  colnames(life_table_male) <- colnames(life_table_female) <- c("age","prob_death","se")
  
  tmp_CA <- rbind(life_table_male[1:111,] %>% 
    mutate(sex="M"), 
    life_table_female[1:111,] %>% 
      mutate(sex="F")) %>% 
    mutate(year=year,
           province = "CA")
  
  tmp_BC <- rbind(life_table_male[1054:1164,] %>% 
                    mutate(sex="M"), 
                  life_table_female[1054:1164,] %>% 
                    mutate(sex="F")) %>% 
    mutate(year=year,
           province = "BC")
  
  CA <- rbind(CA,tmp_CA)
  BC <- rbind(BC,tmp_BC)
}


CA$age <- parse_number(sub("/.*", "", CA$age))
BC$age <- parse_number(sub("/.*", "", BC$age))

write_csv(rbind(CA,BC),"life_table.csv")


# for years > 2020, assume the 2020 life table 