setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

#load("validParameterSets.Rdata")

library(shiny)
library(ggplot2)

#valid_rMechs <- reactionInfo[is.na(reactionInfo$Qvalue) | reactionInfo$Qvalue < 0.1 & !is.na(reactionInfo$Qvalue),]

#unique_rxns = unique(sort(valid_rMechs$reaction))


#valid_rMechsd_rx_info <- rxnf[names(rxnf) %in% valid_rMechs$rMech]

# find all pathways
# point each pathway to all relevent reactions using a list - order as factors with ALL being first
# generate figures and dump out when appropriate


pathwaySet
rxToPW
reactionInfo



ggplotList <- list()

ggplotList$norm <- ggplot(data.frame(x = rnorm(10000, 0, 1)), aes(x = x)) + geom_bar()

ggplotList$pois <- ggplot(data.frame(x = rpois(10000, 5)), aes(x = x)) + geom_bar()


shinyServer(function(input, output) {

  # pathway selected
  pathways_selected <- reactive({input$pathway})
  
  # reactions available
  reactions_available <- reactive({rxToPW$reactionName[rxToPW$pathway == pathwaySet$pathway[pathwaySet$display == pathways_selected()]]})
  output$reactions_available <- renderUI({selectInput("reaction_chosen", "Reaction:", as.list(reactions_available()))})
  
  # which reaction was selected - using common name
  rID = reactive({rxToPW$rID[rxToPW$reactionName == input$reaction_chosen][1]})
  subtypeChoices <- reactive({reactionInfo$Name[reactionInfo$reaction == rID()]})
  output$subforms_available <- renderUI({selectInput("subform_chosen", "Subtype:", as.list(subtypeChoices()))})
  
  #######
  
  
  
  
  # specify which reaction form is desired
  #  
  output$test <- renderText(rID())
  
  output$indicator <- renderText(subtypeChoices())

  

 
  # Generate a plot of the requested variable against mpg and only 
  # include outliers if requested
  output$mpgPlot <- renderPlot({
    boxplot(as.formula(formulaText()), 
            data = mpgData,
            outline = input$outliers)
  })
})
