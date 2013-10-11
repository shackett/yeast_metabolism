library(shiny)


shinyUI(pageWithSidebar(

  headerPanel("Reaction Omics Consistency Analysis (rOCA)"),

  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    
    # Select the pathway of interest (or all reactions)
    selectInput("pathway", "Pathway:", as.list(pathwaySet$display)),
    
    # Select the reaction within this pathway which is to be considered
    uiOutput("reactions_available"),
    
    # Select the reaction subtype of interest
    uiOutput("subforms_available"),
    
    checkboxGroupInput("plots_chosen", "Plots to choose", choices = as.list(names(ggplotList)))
      
    ),
    
    # Partial example
	
    
    #uiOutput("rxn_input")
    
    #		as.list(names(reactionNames)[reactionChoices %in% reaction_selected])),
    
  mainPanel(
    h3(textOutput("indicator")),
    
    h3(textOutput("test")),
    
    plotOutput("P1")
  )
))
