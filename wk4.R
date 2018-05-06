df<-readr::read_csv("daily_SPEC_2014.csv.bz2")


library(dplyr)
df %>% filter(`State Name` == "Wisconsin" & `Parameter Name` == "Bromine PM2.5 LC") %>% 
  summarize(mean = mean(`Arithmetic Mean`, na.rm=TRUE))

df %>% filter(`Parameter Name` %in% c("OC CSN Unadjusted PM2.5 LC TOT",
                                      "EC2 PM2.5 LC",
                                      "Sodium PM2.5 LC",
                                      "Sulfur PM2.5 LC")) %>% group_by(`Parameter Name`) %>% summarize(Mean = mean(`Arithmetic Mean`, na.rm=TRUE))

df %>% filter(`Parameter Name` == "Sulfate PM2.5 LC") %>% 
  filter(`State Code` %in% c("31", "48", "39", "10")) %>%
  filter(`County Code` %in% c("055", "203", "081", "001")) %>%
  filter(`Site Num` %in% c("0019", "0002", "0017", "0003")) %>% 
  group_by(`State Code`, `County Code`, `Site Num`) %>% 
  summarize(Mean = mean(`Arithmetic Mean`, na.rm=TRUE))
           
df %>% filter(`State Name` %in% c("California", "Arizona")) %>% 
  filter(`Parameter Name` == "EC PM2.5 LC TOR") %>% 
  group_by(`State Name`) %>%
  summarize(Mean = mean(`Arithmetic Mean`, na.rm=TRUE))

df %>% filter(Longitude < -100) %>% 
  filter(`Parameter Name` == "OC PM2.5 LC TOR") %>% 
  summarize(Median = median(`Arithmetic Mean`, na.rm=TRUdfE))

# 6
df1 <- readxl::read_excel("aqs_sites.xlsx")

df1 %>% filter(`Location Setting` == "SUBURBAN" & `Land Use` == "RESIDENTIAL") %>% 
  summarize(N=n())

# 7
site <- rename(df1, `Site Num` = `Site Number`) %>%
  select(`State Code`, `County Code`, `Site Num`, `Longitude`, `Land Use`, 
         `Location Setting`)
str(site)
subdata <- mutate(df, `State Code` = as.numeric(`State Code`),
                  `County Code` = as.numeric(`County Code`),
                  `Site Num` = as.numeric(`Site Num`)) %>%
  select(`State Code`, `County Code`, `Site Num`, `Parameter Name`, `Arithmetic Mean`, `Date Local`)

str(subdata) # make sure variables are in the same class
m <- left_join(subdata, site, by = c("State Code", "County Code", "Site Num"))

str(m)
m %>% filter(`Parameter Name` == "EC PM2.5 LC TOR" & `Land Use` == "RESIDENTIAL" &  
               `Location Setting` == "SUBURBAN" & Longitude >= -100) %>%
  group_by("Parameter Name")%>%
  summarize(median=median(`Arithmetic Mean`, na.rm = TRUE))

# 8
subdata %>% left_join(site,by = c("State Code", "County Code", "Site Num")) %>% 
  filter(`Land Use` == "COMMERCIAL",
         `Parameter Name` == "Sulfate PM2.5 LC") %>% 
  mutate(month = lubridate::month(`Date Local`, label=TRUE)) %>% 
  group_by(month) %>%
  summarise(mean = mean(`Arithmetic Mean`, na.rm = TRUE)) %>%
  arrange(desc(mean))

# 9
df %>% filter(`State Code` == "06" & `County Code` == "065" & `Site Num` == "8001" &
              `Parameter Name` %in% c("Sulfate PM2.5 LC", "Total Nitrate PM2.5 LC")) %>% 
  group_by(`Parameter Name`, `Date Local`) %>% 
  select(`State Code`, `County Code`, `Site Num`, `Date Local`, `Parameter Name`,
         `Arithmetic Mean`) %>%
  summarise(mean = mean(`Arithmetic Mean`, na.rm = TRUE)) %>%
  group_by(`Date Local`) %>%
  summarise(Total = sum(mean, na.rm = TRUE)) %>%
  filter(Total > 10)
  
# 10
df %>%
  filter(`Parameter Name` %in% c("Sulfate PM2.5 LC", "Total Nitrate PM2.5 LC")) %>%
  group_by(`State Code`, `County Code`, `Site Num`, `Parameter Name`, `Date Local`) %>%
  select(`State Code`, `County Code`, `Site Num`, `Date Local`, `Parameter Name`,
         `Arithmetic Mean`) %>%
  summarise(mean = mean(`Arithmetic Mean`, na.rm = TRUE)) %>%
  tidyr::spread(`Parameter Name`, mean) %>%
  group_by(`State Code`, `County Code`, `Site Num`) %>%
  summarise(correlation = cor(`Sulfate PM2.5 LC`, `Total Nitrate PM2.5 LC`)) %>%
  arrange(desc(correlation))
