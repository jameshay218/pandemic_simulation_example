#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
#options(shiny.maxRequestSize=1000*1024^2)
shinyUI(


    # Sidebar with a slider input for number of bins
    fluidPage(
        # Application title
        titlePanel("Fitting a deterministic age-structured SEIR model to epi data"),
            tabPanel("Model",
                sidebarPanel(
                    sliderInput("R0",
                                HTML(paste0("Basic reproductive number, R",tags$sub("0"),":")),
                                min = 1,
                                max = 10,
                                value = 2.5,
                                step=0.1),
                    sliderInput("gamma",
                                "Infectious period (days):",
                                min = 1,
                                max = 10,
                                value = 5,
                                step=0.1),
                    sliderInput("alpha",
                                "Latent period (days):",
                                min = 1,
                                max = 10,
                                value = 5,
                                step=0.1),
                   
                    hr(),
                    fluidRow(column(4, h4("Susceptibility")),column(4,h4("Infectivity"))),
                    fluidRow(
                        column(4,sliderInput("age_1_sus","",
                            min = 0,max = 1,value = 0.5,step=0.01,ticks=FALSE)),
                       column(4,sliderInput("age_1_infect","",
                                            min = 0,max = 1,value = 0.25,step=0.01,ticks=FALSE)),
                       column(3, h5("<9"))
                    ),
                    
                    fluidRow(
                        column(4,sliderInput("age_2_sus","",
                                             min = 0,max = 1,value = 0.75,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_2_infect","",
                                             min = 0,max = 1,value = 0.25,step=0.01,ticks=FALSE)),
                        column(3, h5("10-19"))
                    ),
                    fluidRow(
                        column(4,sliderInput("age_3_sus","",
                                             min = 0,max = 1,value = 1,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_3_infect","",
                                             min = 0,max = 1,value = 0.3,step=0.01,ticks=FALSE)),
                        column(3, h5("20-29"))
                    ),
                    fluidRow(
                        column(4,sliderInput("age_4_sus","",
                                             min = 0,max = 1,value = 0.75,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_4_infect","",
                                             min = 0,max = 1,value = 0.45,step=0.01,ticks=FALSE)),
                        column(3, h5("30-39"))
                    ),
                    fluidRow(
                        column(4,sliderInput("age_5_sus","",
                                             min = 0,max = 1,value = 0.6,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_5_infect","",
                                             min = 0,max = 1,value = 0.6,step=0.01,ticks=FALSE)),
                        column(3, h5("40-49"))
                    ),
                    fluidRow(
                        column(4,sliderInput("age_6_sus","",
                                             min = 0,max = 1,value = 0.5,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_6_infect","",
                                             min = 0,max = 1,value = 0.75,step=0.01,ticks=FALSE)),
                        column(3, h5("50-59"))
                    ),
                    fluidRow(
                        column(4,sliderInput("age_7_sus","",
                                             min = 0,max = 1,value = 0.5,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_7_infect","",
                                             min = 0,max = 1,value = 0.8,step=0.01,ticks=FALSE)),
                        column(3, h5("60-69"))
                    ),
                    
                    fluidRow(
                        column(4,sliderInput("age_8_sus","",
                                             min = 0,max = 1,value = 0.5,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_8_infect","",
                                             min = 0,max = 1,value = 0.9,step=0.01,ticks=FALSE)),
                        column(3, h5("70-79"))
                    ),
                    
                    fluidRow(
                        column(4,sliderInput("age_9_sus","",
                                             min = 0,max = 1,value = 0.5,step=0.01,ticks=FALSE)),
                        column(4,sliderInput("age_9_infect","",
                                             min = 0,max = 1,value = 1,step=0.01,ticks=FALSE)),
                        column(3, h5("80+"))
                    ),
                    hr(),
                    fluidRow(
                        column(4,
                               numericInput("N_tot",
                                            "Population size:",
                                            min = 1000,
                                            max = 10000000,
                                            value = 1000000,
                                            step=1000)),
                        column(4,
                               numericInput("symptomatic_proportion",
                                            "Proportion symptomatic:",
                                            min = 0.01,
                                            max = 1,
                                            value = 0.1,
                                            step=0.01))
                    ),
                    fluidRow(
                        column(4,
                               numericInput("reporting_rate",
                                            "Proportion reported:",
                                            min = 0.01,
                                            max = 1,
                                            value = 0.5,
                                            step=0.01)),
                        column(4,
                               numericInput("I0","Seed size (I0):",
                                            min = 1,
                                            max = 10000,
                                            value = 10,
                                            step=1))
                    ),
                    fluidRow(
                        column(4, numericInput("tmin","Min time",0,min=0,max=100)),
                        column(4,numericInput("tmax","Max time",100,min=50,max=500))
                    )
                    
                    ),
                
                
                width=4
                ),
            # Show a plot of the generated distribution
            mainPanel(
                plotOutput("main_plot",width="100%")
                #tabsetPanel(tabPanel("Main",plotOutput("main_plot",width="100%")))
            )
        
    )
)
