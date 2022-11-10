
#path is the path to file to be read
#file that contains beta and allele frequency information
#name is is the string for trait name, i.e "BMI"
#population is the string for population name, i.e. "Non-Finnish European"


#reader will take the given allele frequency and bin them 
#base on different levels of sum in genetic variance and difference in genetic variance



reader_ver5 <- function(path1, name, population){
    
    library(dplyr)
    library(tidyr)
    test  <- read.table(path1, header=T)
    
    #gathering all the necessary columns, and remove any with missing betas
    test <- test %>% select(grep("^AF|BETA", colnames(test))) %>%
                      filter(!is.na(FEMALE_BETA))
    
    
    #Calculate delta var = p*(1-p)*(betafemale ^2 + betamale^2)
    #Calculate delta var quant
    
    #Calculate sum var = p*(1-p)*(betafemale^2 + betamale^2)
    #Calculate sum var quant
    
    mtest5 <- pivot_longer(test, 
                          grep("AF", colnames(test)), 
                          names_to = "POPULATION", 
                          values_to = "AF")
    
    #Treat each population as a factor, separated
    mtest5 <- mtest5 %>% 
              filter(AF <= 0.5) %>%
              filter(AF != 0) %>%
              filter(!is.na(AF)) %>%
              mutate(SUM_VAR = AF*(1-AF)*(FEMALE_BETA^2 + MALE_BETA^2)) %>%
              group_by(POPULATION) %>%
              mutate(SUM_VAR_QUANT = as.factor(ntile(SUM_VAR, 8)))
    
    mtest5 <- mtest5 %>% mutate(DELTA_VAR = AF*(1-AF)*(FEMALE_BETA^2 - MALE_BETA^2)) %>%
              group_by(POPULATION, SUM_VAR_QUANT) %>%
              mutate(DELTA_VAR_QUANT = as.factor(ntile(DELTA_VAR, 30)))
    
    #write.table(mtest5, paste0(name, "_ver5"), row.names=F)
    
    pop <- mtest5 %>% filter(POPULATION == as.character(population))
    
    #Calculate the mean and standard deviation within each of these buckets
    #only doing calcs on minor allele frequencies of nfe population
    
    pop_mean_sd <- pop %>%
         group_by(SUM_VAR_QUANT, DELTA_VAR_QUANT) %>%
         summarize(MEAN_AF = mean(AF, na.rm=TRUE), SD = sd(AF, na.rm=TRUE))
 
    pop_count <- pop %>%
          count(DELTA_VAR_QUANT)
    
    merged_pop <- merge(pop_mean_sd, pop_count) %>% 
              mutate(SE = SD/sqrt(n))
    
    write.table(merged_pop, paste0(name, ".plotdata"))
    return(merged_nfe)

}


#population: name of population as a string
#name: name of trait, as a string
#table: input table generated from above
#Create a plot with mean +/- 2SE of allele frequency for a given trait


plotter_ver5 <- function(table, name, population){ 
    library(ggplot2)

    
    levels(table$SUM_VAR_QUANT) <- list(   "Sum var 0%-12%" = 1,        # Change factor levels
                                           "Sum var 12%-25%" = 2,
                                           "Sum var 25%-37%" = 3,
                                           "Sum var 37%-50%" = 4,
                                           "Sum var 50%-62%" = 5,
                                           "Sum var 62%-75%" = 6,
                                           "Sum var 75%-87%" = 7, 
                                           "Sum var 87%-100%" = 8)
    
    options(repr.plot.width = 18, repr.plot.height=9)
    plot <- table %>%
    ggplot(aes(x=DELTA_VAR_QUANT, y = MEAN_AF, color = SUM_VAR_QUANT, group=SUM_VAR_QUANT)) +
    geom_point(color = 'black', size = 0.8) + 
    geom_line() +
    geom_errorbar(aes(ymin=MEAN_AF - 2*SE, ymax=MEAN_AF + 2*SE), width=.1) + 
    theme_minimal() +
    scale_color_manual( values = c('#A8E1FD', '#A0D8F4', '#85C5E4', 
                                   '#68ACCD', '#5897B6', '#4E81B6', 
                                   '#3B6999', '#2A4E80')) +
    labs(title = paste0("Mean +/- 2SE minor Allele frequency of", population, " population: \n", name)) +
    xlab('Increasing levels of genetic variance difference within sex') + 
    ylab('Minor allele frequency') +
    guides(color = guide_legend(title = "Bins of\nAdditive genetic variance\nwithin sex")) +
    theme(  plot.title = element_text(size = 30),
            plot.subtitle = element_text(size= 18),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.text.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.title = element_text(size = 20,face = "bold"),
            strip.text = element_text(size = 20),
            strip.text.x = element_text(size = 20),
            legend.position='right',
            legend.background = element_rect(color = 'white'),
            legend.text = element_text(size=16),
            legend.title = element_text(size = 20),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.spacing = unit(2, "lines"))

    #Uncomment if need to save the plot or see the plot
    #show(plot)
    #ggsave(paste0(name, "_", population, ".pdf"), width = 18, height = 9)

}
