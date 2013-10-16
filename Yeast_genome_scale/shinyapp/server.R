setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

#load("validParameterSets.Rdata")

library(shiny)
library(ggplot2)
library(gridExtra)

#valid_rMechs <- reactionInfo[is.na(reactionInfo$Qvalue) | reactionInfo$Qvalue < 0.1 & !is.na(reactionInfo$Qvalue),]

#unique_rxns = unique(sort(valid_rMechs$reaction))


#valid_rMechsd_rx_info <- rxnf[names(rxnf) %in% valid_rMechs$rMech]

# find all pathways
# point each pathway to all relevent reactions using a list - order as factors with ALL being first
# generate figures and dump out when appropriate

##### Five files loaded #####
# pathwaySet - number of members in various pathways
# rxToPW - which pathways is a reaction associated with + name
# reactionInfo - reactiion -> rxMechanism + sanitized names
# pathway_plot_list - plots showing information (performance) at a pathway level
# shiny_flux_data - rxMechanism specific plots

load("shinyapp/shinyData.Rdata")


pathwaySet
rxToPW
reactionInfo
pathway_plot_list
shiny_flux_data


ggplotList <- list()

ggplotList$norm <- ggplot(data.frame(x = rnorm(10000, 0, 1)), aes(x = x)) + geom_bar()

ggplotList$pois <- ggplot(data.frame(x = rpois(10000, 5)), aes(x = x)) + geom_bar()


shinyServer(function(input, output) {

  # pathway selected
  pathways_selected <- reactive({input$pathway})
  
  # reactions available
  reactions_available <- reactive({c("PATHWAY INFORMATION", rxToPW$reactionName[rxToPW$pathway == pathwaySet$pathway[pathwaySet$display == pathways_selected()]])})
  output$reactions_available <- renderUI({selectInput("reaction_chosen", "Reaction:", as.list(reactions_available()))})
  
  viewType = reactive({ifelse(input$reaction_chosen == "PATHWAY INFORMATION", "PW", "RX")})
  
  ### If a reaction is chosen ###
  
  rID = reactive({rxToPW$rID[rxToPW$reactionName == input$reaction_chosen][1]})
  subtypeChoices <- reactive({reactionInfo$Name[reactionInfo$reaction == rID()]})
  
  output$subforms_available <- renderUI({selectInput("subform_chosen", "Subtype:", as.list(subtypeChoices()))})
  currentRx <- reactive({reactionInfo$rMech[reactionInfo$reaction == rID()][reactionInfo$Name[reactionInfo$reaction == rID()] == input$subform_chosen]})
  
  RXplots <- reactive({names(shiny_flux_data[[currentRx]])})
  
  ### If a pathway is chosen ###
  
  PWplots <- reactive({names(pathway_plot_list[[pathways_selected()]])})
  
  ### Choose which plots to display ###
  
  plotChoose <- function(view){
    if(view == "PW"){
     return(PWplots()) 
    }else{
     return(RXplots()) 
    }
  }
    
  
  namez <- reactive({plotChoose(viewType())})
  
  output$choosePlot <- renderUI({
    checkboxGroupInput("plots_chosen", "Plots to choose", choices = as.list(displayPlots))
  })
  
  
  
  
  
  plots <- reactive({input$plots_chosen})
  
  plot_subset <- reactive({ggplotList[names(ggplotList) %in% plots()]})
  
  column_number <- reactive({as.numeric(input$col_num)})
  
  # specify which reaction form is desired
    
  output$indicator <- renderText(namez())

  output$test <- renderText(PWplots())
  
  # Generate a plot of the requested variable against mpg and only 
  # include outliers if requested
  
  #grid.arrange(plot1, plot2, ncol=2)
  
  output$P1 <- renderPlot({
    if(length(plot_subset) != 0){
      subList <- plot_subset()
      subList$ncol = column_number()
      
      do.call(grid.arrange,  subList)
    }else{
      return(NULL) 
    }
  })
  
  
})
