setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

load("validParameterSets.Rdata")

library(shiny)

valid_rMechs <- reactionInfo[is.na(reactionInfo$Qvalue) | reactionInfo$Qvalue < 0.1 & !is.na(reactionInfo$Qvalue),]

unique_rxns = unique(sort(valid_rMechs$reaction))


valid_rMechsd_rx_info <- rxnf[names(rxnf) %in% valid_rMechs$rMech]

# find all pathways
# point each pathway to all relevent reactions using a list - order as factors with ALL being first
# generate figures and dump out when appropriate




shinyServer(function(input, output) {

  # which reaction was selected - using common name
  
  reaction_selected <- reactive({input$reaction})
  subrxns <- reactive({valid_rMechs$rMech[valid_rMechs$reaction == reaction_selected()]})
	
  

  # specify which reaction form is desired
  #  
  output$test <- renderText(subrxns())
  
  output$indicator <- renderText({reaction_selected()})

  output$subrxns <- renderUI({selectInput("subrxn", "subrxn", as.list(subrxns()))})
  




 
  # Generate a plot of the requested variable against mpg and only 
  # include outliers if requested
  output$mpgPlot <- renderPlot({
    boxplot(as.formula(formulaText()), 
            data = mpgData,
            outline = input$outliers)
  })
})
