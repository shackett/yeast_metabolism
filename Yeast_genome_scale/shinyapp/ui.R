library(shiny)


shinyUI(pageWithSidebar(

  headerPanel("Reaction Omics Consistency Analysis (rOCA)"),

  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("reaction", "Reaction:",
                as.list(unique_rxns)),
    
    uiOutput("subrxns")            
    ),
    
    # Partial example
	
    
    #uiOutput("rxn_input")
    
    #		as.list(names(reactionNames)[reactionChoices %in% reaction_selected])),
    
  mainPanel(
    h3(textOutput("indicator")),
    
    h3(textOutput("test")),
    
    plotOutput("mpgPlot")
  )
))
