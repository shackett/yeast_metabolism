library(shiny)


shinyUI(pageWithSidebar(

  headerPanel("Reaction Omics Consistency Analysis (rOCA)"),

  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    
    # Select the pathway of interest (or all reactions)
    selectInput("pathway", "Pathway:", as.list(pathwaySet$display)),
    
    uiOutput("pw_check"),
    
    # Allow for saving pathway plots and reaction plots
    textInput('pw_filename', "Pathway Filename"),
    checkboxInput('pw_save', "Save pathway plot"),
    
    
    # Select the reaction within this pathway which is to be considered
    uiOutput("reactions_available"),
    
    # Select the reaction subtype of interest
    uiOutput("subforms_available"),
    
    # Which subset of plots to choose
    uiOutput("rx_check"),
    
    #checkboxGroupInput("plots_chosen", "Plots to choose", choices = as.list(names(ggplotList))),
    
    # How many columns should ggplot figures be displayed in
    selectInput("col_num", "Number of plot columns", c("1", "2", "3", "4"), selected = "2"),
    
    # Allow for saving pathway plots and reaction plots
    textInput('rxn_filename', "Reaction Filename"),
    checkboxInput('rxn_save', "Save reaction plot")
    
    ),
    
  mainPanel(
    plotOutput("PW", height = 1000),
    
    plotOutput("RXinfo", height = 400),
    
    plotOutput("RX", height = 1000)
  )
))
