"Volcano PLot Function

#################

Basilin Benson
University of Washington, basilinb@uw.edu
Copyright (C) 2021 Basilin Benson
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Input parameters:
REQUIRED
    file_name = Thsi is the file generated from extract_pval function
#Input should be a data frame, not file
#For consistency with kimma, suggest df name is 'fdr'

OPTIONAL
  coeff = By default the function will print all coefficients from the linear model. can also specify specific coefficients e.g. c('coeff1', 'coeff2')
  fdr = You can set the FDR correction by default it is 0.2
  n_top_gene = number of up regulated and down regulated genes to label with hgnc symbol and color. default is 10
  mod_gene = you can specify if its modules or genes in you data for the plot. Defaut is gene

#Suggest parameters like
## col_cutoff = fdr cutoff for colors
## label_n = to replace n_top_gene
## label_cutoff = fdr cutoff to label genes (error message if both col_cutoff and label_cutoff are set)

Example
  volcano_funcxn(file_name='results/test.pval_w_hgnc.csv',
             coeff='coeff1',
             fdr='0.05', n_top_gene=15,mod_gene='gene',
             out_dir='./')
"
# All the libraries
volcano_funcxn <- function(file_name, coeff = "all",fdr = 0.2, n_top_gene = 10,mod_gene = "gene", out_dir = "./") {
#Packages do not need to be installed for a function in a custom package
#Similarly in a package, library() is not used. Instead, designate each non-base R function to a given package such as dplyr::select()
  for (package in c('tidyverse', 'ggplot2', 'ggrepel')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }
  }
  options(ggrepel.max.overlaps = Inf)
  # Read in file generated from extract pvalue function for model of interest
  p_val <- read_csv(file_name) %>% arrange(adj.P.Val)
  #create directory
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # If you want to run function for all coefficients
  if(coeff == "all") {
    coeff <- unique(p_val$group)
    coeff <- coeff [! coeff %in% "(Intercept)"]
  }
  else {
  }
    #A lot of code if repeated in this ifelse section
    #Suggest instead using input parameters to set custom column names within the fdr df. Defaults should be those from kimma::fmFit output
  if(mod_gene == "gene"){
    # running for loop for each given coefficient
  for (i in coeff) {
    pval_coef <- subset(p_val, group == i)
    # getting all hgnc symbol for up regulated genes
    up_reg <-  head(subset(pval_coef, FC.group == "up" & adj.P.Val < fdr) %>% arrange(adj.P.Val), n_top_gene) %>% #FC.group may not be present in other user's data. Suggest using the estimate (logFC) to define up/down from pos/neg
              rename(delabel = hgnc_symbol) %>% #Suggest parameter to set what column to pull label info from
              select(geneName, delabel)
    # getting all hgnc symbol for down regulated genes
    down_reg <-  head(subset(pval_coef, FC.group == "down" & adj.P.Val < fdr) %>% arrange(adj.P.Val), n_top_gene) %>%
                rename(delabel = hgnc_symbol) %>%
                select(geneName, delabel)
    # making the dataframe to plot
    pval_plot <- left_join(pval_coef,down_reg, by="geneName") %>% #mutate(ifelse) is a shorter way to add a delabel column to the df
                left_join(.,up_reg,by="geneName") %>%
                mutate(delabel = coalesce(delabel.x,delabel.y)) %>%
                mutate(FC.group = case_when(is.na(delabel) == FALSE ~ FC.group,
                                  is.na(delabel) == TRUE ~ "No")) #No is a confusing label for FC group. These have a fold change, they just aren't colored
    cols <- c("up"="red","down"="blue","No"="black") #Suggest parameter to custom set colors if desired. Defaults maybe should be a little less harsh than these - diffcult to read the text labels
    # running ggplot to make volcano plot
    plt_vol<-ggplot(data=pval_plot, aes(x=logFC, y=-log10(adj.P.Val), col=FC.group, label=delabel)) + 
      #Suggest creating a base ggplot without all aes set then using if statements to add layers like color and labeling only if 
      #they are given as input parameters. This allows the user to use color, label, both, or neither
      geom_point() +
      theme_minimal() +
      geom_text_repel() +
      scale_color_manual(values=cols) +
      geom_vline(xintercept=0, col="black") +
      geom_hline(yintercept= 0, col="black")
      #Suggest vline at fdr cutoff if used
    ggsave(paste(out_dir,i,"_volcano_plot.pdf")) #Save plots to list object within R instead
  }
  }
    else {
      for (i in coeff) {
        pval_coef <- subset(p_val, group == i)
        # getting all hgnc symbol for up regulated genes
        up_reg <-  head(subset(pval_coef, FC.group == "up" & adj.P.Val < fdr) %>% arrange(adj.P.Val), n_top_gene) %>%
          mutate(delabel = geneName) %>%
          select(geneName, delabel)
        # getting all hgnc symbol for down regulated genes
        down_reg <-  head(subset(pval_coef, FC.group == "down" & adj.P.Val < fdr) %>% arrange(adj.P.Val), n_top_gene) %>%
          mutate(delabel = geneName) %>%
          select(geneName, delabel)
        # making the dataframe to plot
        pval_plot <- left_join(pval_coef,down_reg, by="geneName") %>%
          left_join(.,up_reg,by="geneName") %>%
          mutate(delabel = coalesce(delabel.x,delabel.y)) %>%
          mutate(FC.group = case_when(is.na(delabel) == FALSE ~ FC.group,
                                      is.na(delabel) == TRUE ~ "No"))
        cols <- c("up"="red","down"="blue","No"="black")
        # running ggplot to make volcano plot
        plt_vol<-ggplot(data=pval_plot, aes(x=logFC, y=-log10(adj.P.Val), col=FC.group, label=delabel)) +
          geom_point() +
          theme_minimal() +
          geom_text_repel() +
          scale_color_manual(values=cols) +
          geom_vline(xintercept=0, col="black") +
          geom_hline(yintercept= 0, col="black")
        ggsave(paste(out_dir,i,"_volcano_plot.pdf"))
    }
  }
  }

